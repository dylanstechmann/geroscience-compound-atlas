"""Hypothesis-mode generation: harder gates, anti-clone cap, written cards.

Play mode remains `python -m gen.pipeline`. This module does not recommend
ingestion, synthesis, or a personal stack.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from gen.filters import HypothesisGates, evaluate_gates
from gen.ga import MolecularGA
from gen.scorer import CompositeSurrogateScorer

MFPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

DEFAULT_CONFIG = Path("configs/hypothesis.yaml")
NONCLAIM = (
    "Research artifact only. Not a drug, supplement, dose, or protocol. "
    "A high mTOR surrogate score is not an IC50. QED is a historical oral-drug prior."
)


def load_hypothesis_config(path: str | Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def fingerprint(mol: Chem.Mol):
    return MFPGEN.GetFingerprint(mol)


def nearest_tanimoto(mol: Chem.Mol, reference_fps: list) -> float:
    if not reference_fps:
        return 0.0
    fp = fingerprint(mol)
    return float(max(DataStructs.TanimotoSimilarity(fp, ref) for ref in reference_fps))


def training_active_fingerprints(bench_df: pd.DataFrame) -> tuple[list, set[str]]:
    fps = []
    keys: set[str] = set()
    actives = bench_df[bench_df["active"] == 1]
    for _, row in actives.iterrows():
        smi = row.get("canonical_smiles")
        if not isinstance(smi, str):
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fps.append(fingerprint(mol))
        key = row.get("inchikey")
        if isinstance(key, str) and key:
            keys.add(key)
    return fps, keys


def clone_penalty(nearest: float, start: float, weight: float) -> float:
    if nearest <= start:
        return 0.0
    return weight * (nearest - start) / max(1e-6, 1.0 - start)


def write_card(row: pd.Series, out_dir: Path, objective: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{row['card_id']}.md"
    reasons = row.get("reject_reasons") or []
    body = f"""# {row["card_id"]}

{NONCLAIM}

## Objective

{objective.strip()}

## Structure

- SMILES: `{row["smiles"]}`
- InChIKey: `{row["inchikey"]}`

## Surrogate numbers (not measurements)

| Field | Value |
|---|---:|
| frozen mTOR logistic probability | {row["mtor_prob"]:.4f} |
| QED | {row["qed"]:.4f} |
| nearest Tanimoto to a training active | {row["nearest_tanimoto"]:.3f} |
| PAINS flag | {row["has_pains"]} |
| hypothesis reward | {row["reward"]:.4f} |

## How to kill this card

- Measured pChEMBL on CHEMBL2842 below 6 after a real assay.
- Nearest-neighbor Tanimoto was an underestimate of motif copying.
- Property gates passed but the graph is still unsynthesizable or unstable.
- The 1 µM active/inactive cliff on the surrogate does not match the molecule.

## Reject reasons

{", ".join(reasons) if reasons else "none (accepted into archive)"}
"""
    path.write_text(body, encoding="utf-8")
    return path


def run_hypothesis_pipeline(
    config_path: str | Path = DEFAULT_CONFIG,
    benchmark_parquet: str | Path = "artifacts/benchmark_dataset.parquet",
    splits_json: str | Path = "artifacts/splits.json",
    output_dir: str | Path = "hypotheses",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    cfg = load_hypothesis_config(config_path)
    RDLogger.DisableLog("rdApp.*")
    gates = HypothesisGates.from_mapping(cfg.get("gates", {}))
    anti = cfg.get("anti_clone", {})
    scoring = cfg.get("scoring", {})
    ga_cfg = cfg.get("ga", {})
    max_tanimoto = float(anti.get("max_tanimoto_to_training_active", 0.55))
    penalty_start = float(anti.get("penalty_start", 0.40))
    penalty_weight = float(anti.get("penalty_weight", 1.25))
    max_cards = int(ga_cfg.get("max_cards", 25))

    scorer = CompositeSurrogateScorer(
        benchmark_parquet=benchmark_parquet, splits_json=splits_json
    )
    bench_df = pd.read_parquet(benchmark_parquet)
    ref_fps, train_keys = training_active_fingerprints(bench_df)

    actives = bench_df[bench_df["active"] == 1].sort_values(
        by="pchembl_value", ascending=False
    )
    seed_smiles = actives["canonical_smiles"].dropna().head(25).tolist()

    ga = MolecularGA(
        scorer=scorer,
        pop_size=int(ga_cfg.get("pop_size", 40)),
        n_generations=int(ga_cfg.get("n_generations", 20)),
        mutation_rate=float(ga_cfg.get("mutation_rate", 0.7)),
        crossover_rate=float(ga_cfg.get("crossover_rate", 0.3)),
        elite_size=int(ga_cfg.get("elite_size", 5)),
        qed_weight=float(scoring.get("qed_weight", 0.25)),
        seed=int(ga_cfg.get("seed", 42)),
    )
    raw_df, ga_stats = ga.run(seed_smiles=seed_smiles, max_archive_size=400)

    reject_counts: dict[str, int] = {}
    accepted_rows: list[dict[str, Any]] = []

    for _, row in raw_df.iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        nearest = nearest_tanimoto(mol, ref_fps) if mol is not None else 1.0
        ok, reasons = evaluate_gates(mol, gates, has_pains=bool(row.get("has_pains")))
        if nearest > max_tanimoto:
            ok = False
            reasons = list(reasons) + ["too_close_to_training_active"]
        try:
            ikey = str(row["inchikey"])
        except Exception:
            ikey = ""
        if ikey in train_keys:
            ok = False
            reasons = list(reasons) + ["seen_in_training"]

        for reason in reasons:
            reject_counts[reason] = reject_counts.get(reason, 0) + 1

        reward = float(row["reward"]) - clone_penalty(nearest, penalty_start, penalty_weight)
        record = {
            "smiles": row["smiles"],
            "inchikey": ikey,
            "generation": int(row.get("generation", -1)),
            "mtor_prob": float(row["mtor_prob"]),
            "qed": float(row["qed"]),
            "has_pains": bool(row["has_pains"]),
            "nearest_tanimoto": nearest,
            "reward": max(0.0, reward),
            "accepted": ok,
            "reject_reasons": ";".join(reasons),
        }
        if ok:
            accepted_rows.append(record)

    accepted = pd.DataFrame(accepted_rows)
    if not accepted.empty:
        accepted = accepted.sort_values(by="reward", ascending=False).drop_duplicates(
            subset=["inchikey"]
        )
        accepted = accepted.head(max_cards).reset_index(drop=True)
        accepted.insert(0, "card_id", [f"HYP-{i+1:03d}" for i in range(len(accepted))])
    else:
        accepted = pd.DataFrame(
            columns=[
                "card_id",
                "smiles",
                "inchikey",
                "generation",
                "mtor_prob",
                "qed",
                "has_pains",
                "nearest_tanimoto",
                "reward",
                "accepted",
                "reject_reasons",
            ]
        )

    out = Path(output_dir)
    stamp = date.today().isoformat()
    run_dir = out / f"{stamp}-mtor-hypothesis"
    cards_dir = run_dir / "cards"
    run_dir.mkdir(parents=True, exist_ok=True)
    accepted.to_csv(run_dir / "accepted.csv", index=False)

    readme = run_dir / "README.md"
    readme.write_text(
        "\n".join(
            [
                f"# Hypothesis run {stamp}",
                "",
                NONCLAIM,
                "",
                cfg.get("objective", "").strip(),
                "",
                f"Accepted cards: {len(accepted)}",
                f"GA archive before gates: {len(raw_df)}",
                f"Anti-clone Tanimoto cap: {max_tanimoto}",
                "",
                "Reject counts:",
                *[f"- {k}: {v}" for k, v in sorted(reject_counts.items(), key=lambda kv: -kv[1])],
                "",
            ]
        ),
        encoding="utf-8",
    )

    for _, row in accepted.iterrows():
        write_card(row, cards_dir, cfg.get("objective", ""))

    metrics = {
        "mode": "hypothesis",
        "date": stamp,
        "ga_archive_size": int(len(raw_df)),
        "accepted_cards": int(len(accepted)),
        "max_tanimoto_to_training_active": max_tanimoto,
        "reject_counts": reject_counts,
        "ga_stats": ga_stats,
        "nonclaim": NONCLAIM,
    }
    (run_dir / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    return accepted, metrics


def main() -> None:
    accepted, metrics = run_hypothesis_pipeline()
    print(json.dumps({"accepted_cards": metrics["accepted_cards"], "rejects": metrics["reject_counts"]}, indent=2))
    print(f"wrote {len(accepted)} cards")


if __name__ == "__main__":
    main()
