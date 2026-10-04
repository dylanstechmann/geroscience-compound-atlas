"""Prespecified ChEMBL source-document/structure holdout with mixed assay contexts, not an aging assay."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urljoin, urlparse

import httpx
import numpy as np
import pandas as pd
import sklearn
from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from atlas.features import compute_morgan_fingerprint, compute_rdkit_descriptors
from atlas.normalize import canonicalize_smiles
from bench.config import load_bench_config
from bench.data import activity_labels
from bench.models import prepare_feature_matrices
from bench.neighborhood import maximum_training_tanimoto
from bench.neighborhood_report import _validate_scaffold_partition

API = "https://www.ebi.ac.uk/chembl/api/data/"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_collection(client, resource: str, key: str, params: dict, raw_dir: Path):
    """Require complete, unique, stable-count pagination and preserve every response."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    url = str(httpx.URL(API + resource + ".json", params={**params, "limit": 1000}))
    rows, receipts, seen = [], [], set()
    expected = None
    while url:
        if url in seen or urlparse(url).hostname != "www.ebi.ac.uk":
            raise ValueError("invalid or repeated ChEMBL pagination URL")
        seen.add(url)
        response = client.get(url)
        response.raise_for_status()
        body = response.content
        page = json.loads(body)
        meta = page["page_meta"]
        total = meta["total_count"]
        if isinstance(total, bool) or not isinstance(total, int) or total < 0:
            raise ValueError("invalid ChEMBL total count")
        if expected is None:
            expected = total
        if total != expected or meta["offset"] != len(rows):
            raise ValueError("ChEMBL pagination count changed or offset mismatched")
        batch = page[key]
        if not isinstance(batch, list) or (not batch and len(rows) != expected):
            raise ValueError("incomplete ChEMBL page")
        name = f"{resource}_{len(receipts):04d}.json"
        destination = raw_dir / name
        destination.write_bytes(body)
        receipts.append(
            {
                "path": name,
                "url": url,
                "sha256": sha256(destination),
                "bytes": len(body),
                "rows": len(batch),
            }
        )
        rows.extend(batch)
        next_url = meta["next"]
        url = urljoin("https://www.ebi.ac.uk", next_url) if next_url else None
        if len(rows) > expected or len(receipts) > 100:
            raise ValueError("unexpected or unbounded ChEMBL pagination")
    identity = {"activity": "activity_id", "assay": "assay_chembl_id"}[resource]
    if len(rows) != expected or len({row[identity] for row in rows}) != len(rows):
        raise ValueError("incomplete or duplicate ChEMBL records")
    return rows, receipts


def fetch_snapshot(raw_dir: Path) -> dict:
    if (raw_dir / "receipt.json").exists():
        raise ValueError("snapshot exists; use it unchanged or choose a new directory")
    with httpx.Client(timeout=90, follow_redirects=True) as client:
        status = client.get(API + "status.json")
        status.raise_for_status()
        raw_dir.mkdir(parents=True, exist_ok=True)
        (raw_dir / "status.json").write_bytes(status.content)
        activities, activity_pages = fetch_collection(
            client,
            "activity",
            "activities",
            {
                "target_chembl_id": "CHEMBL2842",
                "standard_type": "IC50",
                "pchembl_value__isnull": "false",
                "order_by": "activity_id",
            },
            raw_dir,
        )
        assays, assay_pages = fetch_collection(
            client,
            "assay",
            "assays",
            {"target_chembl_id": "CHEMBL2842", "order_by": "assay_chembl_id"},
            raw_dir,
        )
    receipt = {
        "retrieved_utc": datetime.now(UTC).isoformat(),
        "database": status.json(),
        "activity_count": len(activities),
        "assay_count": len(assays),
        "status_sha256": sha256(raw_dir / "status.json"),
        "pages": activity_pages + assay_pages,
    }
    (raw_dir / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def load_snapshot(raw_dir: Path):
    receipt = json.loads((raw_dir / "receipt.json").read_bytes())
    if sha256(raw_dir / "status.json") != receipt["status_sha256"]:
        raise ValueError("source status hash mismatch")
    collections = {"activity": [], "assay": []}
    for entry in receipt["pages"]:
        path = raw_dir / entry["path"]
        if path.resolve().parent != raw_dir.resolve():
            raise ValueError("snapshot path must be a direct child")
        body = path.read_bytes()
        if len(body) != entry["bytes"] or hashlib.sha256(body).hexdigest() != entry["sha256"]:
            raise ValueError("source response hash mismatch")
        resource = entry["path"].split("_")[0]
        page = json.loads(body)
        collections[resource].extend(page["activities" if resource == "activity" else "assays"])
    for resource, id_key in (("activity", "activity_id"), ("assay", "assay_chembl_id")):
        rows = collections[resource]
        if len(rows) != receipt[resource + "_count"] or len({r[id_key] for r in rows}) != len(rows):
            raise ValueError("snapshot count or record identities mismatch")
    return collections["activity"], collections["assay"], receipt


def exclusion_reason(activity: dict, assay: dict | None) -> str | None:
    if assay is None:
        return "missing_assay"
    checks = {
        "target_chembl_id": "CHEMBL2842",
        "assay_type": "B",
        "confidence_score": 9,
        "relationship_type": "D",
        "assay_organism": "Homo sapiens",
    }
    for key, value in checks.items():
        if assay.get(key) != value:
            return "assay_" + key
    if assay.get("variant_sequence") is not None or activity.get("assay_variant_mutation"):
        return "variant"
    description = str(assay.get("description", "")).lower()
    description = description.replace("fkbp12-independent", "").replace("fkbp12 independent", "")
    if any(term in description for term in ("fkbp", "fk506 binding", "proliferation", "viability")):
        return "assay_description_scope"
    for key, value in {
        "target_chembl_id": "CHEMBL2842",
        "assay_type": "B",
        "standard_type": "IC50",
        "standard_relation": "=",
        "standard_units": "nM",
        "standard_flag": 1,
        "potential_duplicate": 0,
    }.items():
        if activity.get(key) != value:
            return "activity_" + key
    if activity.get("data_validity_comment") not in (None, "Manually validated"):
        return "data_validity"
    if not activity.get("document_chembl_id") or activity.get("document_chembl_id") != assay.get(
        "document_chembl_id"
    ):
        return "document_mismatch"
    try:
        pvalue, nm = float(activity["pchembl_value"]), float(activity["standard_value"])
        if (
            not np.isfinite([pvalue, nm]).all()
            or nm <= 0
            or abs(pvalue - (9 - np.log10(nm))) > 0.03
        ):
            return "invalid_potency"
    except (KeyError, TypeError, ValueError):
        return "invalid_potency"
    return None


def qualify(activities, assays, original_records, original_frame):
    lookup = {a["assay_chembl_id"]: a for a in assays}
    documents = {r.get("document_chembl_id") for r in original_records}
    molecule_ids = {
        r.get(k)
        for r in original_records
        for k in ("molecule_chembl_id", "parent_molecule_chembl_id")
    }
    molecule_ids.discard(None)
    connectivity = {key[:14] for key in original_frame["inchikey"]}
    counts, accepted = Counter(), []
    for row in activities:
        reason = exclusion_reason(row, lookup.get(row.get("assay_chembl_id")))
        if reason is None and row.get("document_chembl_id") in documents:
            reason = "original_document"
        if reason is None and any(
            row.get(k) in molecule_ids for k in ("molecule_chembl_id", "parent_molecule_chembl_id")
        ):
            reason = "original_molecule"
        if reason is None:
            try:
                smiles, key = canonicalize_smiles(row.get("canonical_smiles"))
            except (ValueError, TypeError):
                reason = "invalid_structure"
            else:
                if key[:14] in connectivity:
                    reason = "original_connectivity"
        if reason:
            counts[reason] += 1
            continue
        accepted.append(
            {
                "connectivity": key[:14],
                "inchikey": key,
                "canonical_smiles": smiles,
                "molecule_chembl_id": row["molecule_chembl_id"],
                "document_chembl_id": row["document_chembl_id"],
                "assay_chembl_id": row["assay_chembl_id"],
                "activity_id": row["activity_id"],
                "pchembl_value": float(row["pchembl_value"]),
            }
        )
    if not accepted:
        return pd.DataFrame(), dict(counts)
    accepted_frame = pd.DataFrame(accepted)
    rows = []
    for _, group in accepted_frame.groupby(["connectivity", "assay_chembl_id"], sort=True):
        row = group.sort_values("activity_id").iloc[0].to_dict()
        row["pchembl_value"] = float(group["pchembl_value"].median())
        row["activity_ids"] = sorted(group["activity_id"].tolist())
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        row["murcko_scaffold"] = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        row.update(compute_rdkit_descriptors(mol))
        rows.append(row)
    counts["accepted_activity_rows"] = len(accepted)
    counts["accepted_molecule_assay_rows"] = len(rows)
    return pd.DataFrame(rows), dict(counts)


def metrics(labels, probabilities):
    both = len(np.unique(labels)) == 2
    return {
        "n": len(labels),
        "positive": int(np.sum(labels)),
        "auroc": float(roc_auc_score(labels, probabilities)) if both else None,
        "average_precision": float(average_precision_score(labels, probabilities))
        if both
        else None,
        "brier": float(brier_score_loss(labels, probabilities)),
    }


def document_bootstrap(labels, probabilities, baseline, documents):
    unique = np.unique(documents)
    if len(unique) < 2:
        return None
    errors = (probabilities - labels) ** 2 - (baseline - labels) ** 2
    # Resample whole documents; primary pooled estimand weights measured rows.
    sums = np.array([errors[documents == d].sum() for d in unique])
    counts = np.array([np.sum(documents == d) for d in unique])
    rng = np.random.default_rng(20261004)
    indexes = rng.integers(0, len(unique), size=(1000, len(unique)))
    pooled = sums[indexes].sum(axis=1) / counts[indexes].sum(axis=1)
    macro = (sums / counts)[indexes].mean(axis=1)
    return {
        "unit": "document",
        "n_documents": len(unique),
        "estimand": "pooled measured-row Brier difference model minus prevalence",
        "difference": float(errors.mean()),
        "percentile_95": np.quantile(pooled, [0.025, 0.975]).tolist(),
        "secondary_equal_document": {
            "difference": float((sums / counts).mean()),
            "percentile_95": np.quantile(macro, [0.025, 0.975]).tolist(),
        },
        "limitation": "documents can share chemistry/labs; intervals do not establish independent biology",
    }


def evaluate(repo: Path, raw_dir: Path, study_dir: Path):
    plan = json.loads((study_dir / "plan.json").read_bytes())
    for path, digest in plan["frozen_inputs"].items():
        if sha256(repo / path) != digest:
            raise ValueError("frozen input changed: " + path)
    config = load_bench_config(repo / "configs/bench.yaml")
    original = pd.read_parquet(repo / "artifacts/benchmark_dataset.parquet")
    records = json.loads((repo / "data/interim/chembl_bench_raw.json").read_bytes())
    activities, assays, receipt = load_snapshot(raw_dir)
    external, counts = qualify(activities, assays, records, original)
    report = {
        "plan_sha256": sha256(study_dir / "plan.json"),
        "source_receipt": receipt,
        "qualification": counts,
        "original_training_assay_qualification": dict(
            Counter(
                exclusion_reason(
                    r, {a["assay_chembl_id"]: a for a in assays}.get(r["assay_chembl_id"])
                )
                or "qualified"
                for r in records
            )
        ),
        "limitations": plan["limitations"],
        "software": {
            "numpy": np.__version__,
            "sklearn": sklearn.__version__,
            "rdkit": rdBase.rdkitVersion,
        },
        "results": {},
    }
    if not external.empty:
        assay_lookup = {a["assay_chembl_id"]: a for a in assays}
        format_labels = [
            assay_lookup[a].get("bao_label") or "unspecified" for a in external["assay_chembl_id"]
        ]
        report["accepted_assay_formats"] = {
            "molecule_assay_rows": dict(Counter(format_labels)),
            "unique_assays": dict(
                Counter(
                    assay_lookup[a].get("bao_label") or "unspecified"
                    for a in set(external["assay_chembl_id"])
                )
            ),
        }
        report["accepted_assay_metadata"] = [
            {
                k: assay_lookup[a].get(k)
                for k in (
                    "assay_chembl_id",
                    "document_chembl_id",
                    "description",
                    "assay_type",
                    "confidence_score",
                    "relationship_type",
                    "assay_organism",
                    "bao_format",
                    "bao_label",
                )
            }
            for a in sorted(set(external["assay_chembl_id"]))
        ]
    report["reporting_corrections"] = [
        "ChEMBL B/confidence9/D metadata qualify target assignment, not exclusively biochemical assay measurement; accepted assays include cellular pathway phosphorylation. The frozen plan scope label was too narrow; eligibility and measurements were not changed.",
        "The primary uncertainty analysis follows the frozen pooled-Brier document-cluster plan. Equal-document Brier is retained only as a secondary estimand.",
        "Independent review identified a false FKBP12-independent keyword rejection. The implementation was corrected and the original cohort/results preserved in pre_correction_results.json. Other FKBP mentions remain conservatively excluded; source eligibility is not adjusted against measured labels.",
    ]
    report["limitations"].append(
        "Accepted target-assigned IC50 assays mix biochemical and cellular pathway contexts; none establish senolysis or rejuvenation."
    )
    if external.empty:
        report["status"] = "no_eligible_external_data"
    else:
        report["status"] = "external_target_assigned_mtor_assay_evaluation"
        original["active"] = activity_labels(
            original["pchembl_value"], config.thresholds.pchembl_active
        )
        external["active"] = activity_labels(
            external["pchembl_value"], config.thresholds.pchembl_active
        )
        frame = pd.concat([original, external], ignore_index=True)
        frame["fingerprint"] = [
            compute_morgan_fingerprint(
                Chem.MolFromSmiles(s), radius=config.features.radius, n_bits=config.features.n_bits
            )
            for s in frame["canonical_smiles"]
        ]
        splits = json.loads((repo / "artifacts/splits.json").read_bytes())
        test = np.arange(len(original), len(frame))
        documents = external["document_chembl_id"].to_numpy()
        y = external["active"].to_numpy()
        predictions = external[
            ["connectivity", "document_chembl_id", "assay_chembl_id", "active"]
        ].copy()
        for seed in config.split.seeds:
            parts = _validate_scaffold_partition(original, splits["scaffold"][str(seed)])
            X, train_y, _, _, external_X, _ = prepare_feature_matrices(
                frame,
                parts["train"],
                parts["val"],
                test,
                include_descriptors=config.features.include_descriptors,
            )
            baseline = np.full(len(y), train_y.mean())
            similarity = maximum_training_tanimoto(
                np.asarray(frame["fingerprint"].tolist()), parts["train"], test
            )
            params = config.model.hgb_params
            models = {
                "logistic": LogisticRegression(
                    C=1.0, max_iter=1000, random_state=seed, solver="lbfgs"
                ),
                "hgb": HistGradientBoostingClassifier(
                    max_iter=params.max_iter,
                    learning_rate=params.learning_rate,
                    min_samples_leaf=params.min_samples_leaf,
                    random_state=seed,
                ),
            }
            result = {
                "training_prevalence_baseline": metrics(y, baseline),
                "similarity": {
                    "median": float(np.median(similarity)),
                    "max": float(similarity.max()),
                    "n_at_most_0_5": int((similarity <= 0.5).sum()),
                },
                "models": {},
            }
            predictions[f"similarity_seed_{seed}"] = similarity
            for name, model in models.items():
                model.fit(X, train_y)
                probabilities = model.predict_proba(external_X)[:, 1]
                predictions[f"{name}_seed_{seed}"] = probabilities
                assays_ids = external["assay_chembl_id"].to_numpy()
                result["models"][name] = {
                    **metrics(y, probabilities),
                    "assay_macro_brier": float(
                        np.mean(
                            [
                                brier_score_loss(y[assays_ids == a], probabilities[assays_ids == a])
                                for a in np.unique(assays_ids)
                            ]
                        )
                    ),
                    "document_bootstrap": document_bootstrap(y, probabilities, baseline, documents),
                    "by_document": {
                        d: metrics(y[documents == d], probabilities[documents == d])
                        for d in np.unique(documents)
                    },
                }
            report["results"][str(seed)] = result
        external.drop(columns=["activity_id"]).to_json(
            study_dir / "qualified_rows.json", orient="records", indent=2
        )
        predictions.to_csv(study_dir / "predictions.csv", index=False)
    (study_dir / "results.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["fetch", "evaluate"])
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw/chembl_external_2026-10-04"))
    parser.add_argument(
        "--study-dir", type=Path, default=Path("studies/chembl_external_2026-10-04")
    )
    args = parser.parse_args()
    if args.command == "fetch":
        receipt = fetch_snapshot(args.raw_dir)
        print(
            json.dumps({k: receipt[k] for k in ("activity_count", "assay_count", "retrieved_utc")})
        )
    else:
        result = evaluate(args.repo, args.raw_dir, args.study_dir)
        print(json.dumps({"status": result["status"], "qualification": result["qualification"]}))


if __name__ == "__main__":
    main()
