"""Benchmark dataset curation: mTOR longevity kinase bioactivities from ChEMBL."""

import json
import logging
from pathlib import Path
from typing import Any

import httpx
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from atlas.features import compute_morgan_fingerprint, compute_rdkit_descriptors
from atlas.normalize import canonicalize_smiles

logger = logging.getLogger(__name__)

CHEMBL_ACTIVITY_URL = "https://www.ebi.ac.uk/chembl/api/data/activity"


def validate_activity_threshold(threshold: float) -> float:
    """Reject invalid thresholds before reading caches or contacting ChEMBL."""
    try:
        value = float(threshold)
    except (TypeError, ValueError) as exc:
        raise ValueError("pchembl_threshold must be a finite number") from exc
    if isinstance(threshold, bool) or not np.isfinite(value):
        raise ValueError("pchembl_threshold must be a finite number")
    return value


def activity_labels(values: pd.Series, threshold: float) -> pd.Series:
    """Derive labels from the stored measurements, never from cached labels."""
    threshold = validate_activity_threshold(threshold)
    try:
        numeric = pd.to_numeric(values, errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError("Benchmark pchembl_value measurements must be finite numbers") from exc
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError("Benchmark pchembl_value measurements must be finite numbers")
    return (numeric >= threshold).astype(int)


def save_benchmark_tables(
    df: pd.DataFrame, output_parquet: str | Path, *, save_fingerprints: bool = False
) -> None:
    """Save cache artifacts before exporting CSV; cached relabeling reuses fingerprints."""
    out_file = Path(output_parquet)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    table = df.drop(columns=["fingerprint"])
    table.to_parquet(out_file, index=False)
    if save_fingerprints:
        # A failed CSV export must leave the new rows paired with their new fingerprints.
        fp_array = np.array(df["fingerprint"].tolist(), dtype=np.float32)
        fp_path = out_file.with_suffix(".fingerprints.npy")
        np.save(fp_path, fp_array)
        logger.info("Saved fingerprint matrix (%s) to %s", fp_array.shape, fp_path)
    csv_file = Path("data/processed/benchmark_dataset.csv")
    csv_file.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(csv_file, index=False)


def fetch_chembl_mtor_records(
    cache_path: str | Path = "data/interim/chembl_bench_raw.json",
    target_chembl_id: str = "CHEMBL2842",
    target_count: int = 600,
) -> list[dict[str, Any]]:
    """Fetch structured IC50 bioactivity records for mTOR from ChEMBL with local caching.

    Args:
        cache_path: Path to disk cache file.
        target_chembl_id: ChEMBL target identifier (mTOR = CHEMBL2842).
        target_count: Target number of records to retrieve.

    Returns:
        List of raw activity dictionaries.
    """
    cache_file = Path(cache_path)
    if cache_file.exists():
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                cached_data = json.load(f)
                if len(cached_data) >= target_count:
                    logger.info(
                        "Loaded %d cached ChEMBL mTOR records from %s", len(cached_data), cache_file
                    )
                    return cached_data
        except (json.JSONDecodeError, OSError) as exc:
            logger.warning("Failed to load cache from %s: %s", cache_file, exc)

    logger.info("Querying ChEMBL API for mTOR (%s) IC50 activities...", target_chembl_id)
    headers = {"User-Agent": "GeroscienceAtlas/0.1.0 (academic research)"}
    records: list[dict[str, Any]] = []
    limit = 100
    offset = 0

    while len(records) < target_count:
        url = (
            f"{CHEMBL_ACTIVITY_URL}?target_chembl_id={target_chembl_id}"
            f"&standard_type=IC50&pchembl_value__isnull=false&limit={limit}&offset={offset}&format=json"
        )
        try:
            with httpx.Client(timeout=20.0) as client:
                res = client.get(url, headers=headers)
            if res.status_code != 200:
                logger.warning("ChEMBL returned status %d at offset %d", res.status_code, offset)
                break
            data = res.json()
            activities = data.get("activities", [])
            if not activities:
                break
            records.extend(activities)
            offset += limit
            logger.info("Retrieved %d/%d records...", len(records), target_count)
        except (httpx.HTTPError, OSError) as exc:
            logger.warning("Error fetching ChEMBL activities at offset %d: %s", offset, exc)
            break

    cache_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(records, f)
    except (TypeError, OSError) as exc:
        logger.warning("Failed to save cache to %s: %s", cache_file, exc)

    return records


def prepare_benchmark_dataset(
    cache_path: str | Path = "data/interim/chembl_bench_raw.json",
    output_parquet: str | Path = "artifacts/benchmark_dataset.parquet",
    pchembl_threshold: float = 6.0,
    target_count: int = 600,
) -> pd.DataFrame:
    """Process raw ChEMBL activities into a normalized, scaffold-annotated benchmark dataset.

    Args:
        cache_path: Path to cached raw records.
        output_parquet: Destination path for frozen dataset parquet file.
        pchembl_threshold: Inclusive cutoff for the active label.
        target_count: Number of records to curate.

    Returns:
        Curated benchmark DataFrame.
    """
    pchembl_threshold = validate_activity_threshold(pchembl_threshold)
    raw_records = fetch_chembl_mtor_records(cache_path=cache_path, target_count=target_count)
    clean_rows = []

    for rec in raw_records:
        raw_smiles = rec.get("canonical_smiles")
        pchembl_str = rec.get("pchembl_value")
        mol_id = rec.get("molecule_chembl_id")

        if not raw_smiles or not pchembl_str:
            continue

        try:
            pchembl = float(pchembl_str)
        except (ValueError, TypeError):
            continue

        try:
            can_smiles, inchikey = canonicalize_smiles(raw_smiles)
        except (ValueError, TypeError):
            continue

        # Extract Bemis-Murcko scaffold
        mol = Chem.MolFromSmiles(can_smiles)
        if mol is None:
            continue

        try:
            scaffold_smiles = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        except (ValueError, RuntimeError):
            scaffold_smiles = ""

        # Compute descriptors
        descriptors = compute_rdkit_descriptors(mol)
        fp = compute_morgan_fingerprint(mol, radius=2, n_bits=2048)

        clean_rows.append(
            {
                "molecule_chembl_id": mol_id,
                "inchikey": inchikey,
                "canonical_smiles": can_smiles,
                "murcko_scaffold": scaffold_smiles,
                "pchembl_value": pchembl,
                "active": int(pchembl >= pchembl_threshold),
                "fingerprint": fp,
                **descriptors,
            }
        )

    df = pd.DataFrame(clean_rows)

    # Deduplicate by InChIKey (mean pchembl_value across replicates)
    dedup_rows = []
    for inchikey, group in df.groupby("inchikey"):
        mean_pchembl = float(group["pchembl_value"].mean())
        first_row = group.iloc[0].to_dict()
        # Preserve precision so a reload cannot move a measurement across the cutoff.
        first_row["pchembl_value"] = mean_pchembl
        dedup_rows.append(first_row)

    final_df = pd.DataFrame(dedup_rows)
    final_df["active"] = activity_labels(final_df["pchembl_value"], pchembl_threshold)

    save_benchmark_tables(final_df, output_parquet, save_fingerprints=True)

    logger.info(
        "Benchmark dataset ready: %d unique molecules (%d active, %d inactive). Scaffolds: %d unique.",
        len(final_df),
        final_df["active"].sum(),
        len(final_df) - final_df["active"].sum(),
        final_df["murcko_scaffold"].nunique(),
    )

    return final_df
