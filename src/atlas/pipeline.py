"""Phase 1 Pipeline: Ingest watchlist, resolve to PubChem identifiers, featurize, and export."""

import argparse
import logging
from pathlib import Path

import pandas as pd

from atlas.features import compute_rdkit_descriptors
from atlas.models import Compound, Modality
from atlas.normalize import canonicalize_smiles
from atlas.pubchem import PubChemClient
from viz.coverage import plot_resolution_coverage

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("atlas.pipeline")


def run_phase1_pipeline(
    watchlist_path: str | Path = "data/watchlist/seed_compounds.csv",
    cache_path: str | Path = "data/interim/pubchem_cache.json",
    parquet_output: str | Path = "artifacts/compounds.parquet",
    figure_output: str | Path = "figures/resolution_coverage.png",
) -> pd.DataFrame:
    """Execute Phase 1 resolution and featurization pipeline.

    Args:
        watchlist_path: Path to input watchlist CSV.
        cache_path: Path to PubChem cache JSON.
        parquet_output: Destination path for frozen compounds parquet file.
        figure_output: Destination path for coverage diagnostic figure.

    Returns:
        Processed DataFrame containing resolved records and descriptors.
    """
    watchlist_file = Path(watchlist_path)
    if not watchlist_file.exists():
        raise FileNotFoundError(f"Watchlist file not found: {watchlist_file}")

    logger.info("Reading watchlist from %s...", watchlist_file)
    watchlist_df = pd.read_csv(watchlist_file)

    client = PubChemClient(cache_path=cache_path)
    records = []

    for _, row in watchlist_df.iterrows():
        raw_name = str(row["raw_name"]).strip()
        raw_modality = str(row.get("modality", "small_molecule")).strip().lower()
        notes = str(row.get("notes", "")).strip()

        try:
            modality = Modality(raw_modality)
        except ValueError:
            modality = Modality.other

        logger.info("Resolving '%s' (%s)...", raw_name, modality.value)
        pubchem_res = client.resolve_name(raw_name)

        if pubchem_res and pubchem_res.get("resolved"):
            cid = pubchem_res.get("cid")
            raw_smiles = pubchem_res.get("canonical_smiles")
            pubchem_key = pubchem_res.get("inchikey")

            # Canonicalize and double-verify InChIKey via RDKit
            try:
                can_smiles, computed_key = canonicalize_smiles(raw_smiles)
                inchikey = computed_key or pubchem_key

                # Validate Pydantic model
                Compound(
                    name=pubchem_res.get("title") or raw_name,
                    synonyms=[raw_name],
                    cid=cid,
                    inchikey=inchikey,
                    canonical_smiles=can_smiles,
                    modality=modality,
                )

                # Compute RDKit descriptors
                descriptors = compute_rdkit_descriptors(can_smiles)
                record = {
                    "raw_name": raw_name,
                    "resolved_name": pubchem_res.get("title") or raw_name,
                    "cid": cid,
                    "inchikey": inchikey,
                    "canonical_smiles": can_smiles,
                    "modality": modality.value,
                    "resolved": True,
                    "notes": notes,
                    **descriptors,
                }
            except (ValueError, TypeError, RuntimeError) as exc:
                logger.warning("RDKit normalization failed for '%s': %s", raw_name, exc)
                record = {
                    "raw_name": raw_name,
                    "resolved_name": pubchem_res.get("title") or raw_name,
                    "cid": cid,
                    "inchikey": pubchem_key,
                    "canonical_smiles": raw_smiles,
                    "modality": modality.value,
                    "resolved": True,
                    "notes": f"RDKit parse warning: {exc}; {notes}".strip(),
                    **compute_rdkit_descriptors(None),
                }
        else:
            logger.warning("Failed to resolve '%s'", raw_name)
            record = {
                "raw_name": raw_name,
                "resolved_name": None,
                "cid": None,
                "inchikey": None,
                "canonical_smiles": None,
                "modality": modality.value,
                "resolved": False,
                "notes": notes,
                **compute_rdkit_descriptors(None),
            }

        records.append(record)

    df_out = pd.DataFrame(records)

    # Save outputs
    out_parquet = Path(parquet_output)
    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_parquet(out_parquet, index=False)
    logger.info("Saved %d compound records to %s", len(df_out), out_parquet)

    # Also save CSV copy in data/processed
    csv_copy = Path("data/processed/compounds.csv")
    csv_copy.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(csv_copy, index=False)

    # Generate coverage plot
    fig_path = plot_resolution_coverage(df_out, output_path=figure_output)
    logger.info("Saved coverage plot to %s", fig_path)

    # Log summary
    resolved_count = df_out["resolved"].sum()
    total_count = len(df_out)
    logger.info(
        "Phase 1 Complete: %d/%d (%.1f%%) compounds successfully resolved.",
        resolved_count,
        total_count,
        (resolved_count / total_count * 100) if total_count > 0 else 0,
    )

    return df_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Phase 1 Resolution Pipeline")
    parser.add_argument("--watchlist", default="data/watchlist/seed_compounds.csv")
    parser.add_argument("--cache", default="data/interim/pubchem_cache.json")
    parser.add_argument("--output", default="artifacts/compounds.parquet")
    parser.add_argument("--figure", default="figures/resolution_coverage.png")
    args = parser.parse_args()

    run_phase1_pipeline(
        watchlist_path=args.watchlist,
        cache_path=args.cache,
        parquet_output=args.output,
        figure_output=args.figure,
    )
