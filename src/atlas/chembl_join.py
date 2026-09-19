"""Phase 2 Pipeline: Join compounds to ChEMBL bioactivities and build curated evidence graph."""

import argparse
import logging
from pathlib import Path
from typing import Any

import pandas as pd

from atlas.chembl import ChEMBLClient
from atlas.graph import EvidenceGraph
from atlas.models import EvidenceEdge, EvidenceGrade, HallmarkSlug
from viz.hallmarks import plot_chembl_and_hallmark_coverage

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("atlas.chembl_join")


def run_phase2_pipeline(
    compounds_parquet: str | Path = "artifacts/compounds.parquet",
    curated_edges_csv: str | Path = "data/curated/curated_evidence.csv",
    cache_path: str | Path = "data/interim/chembl_cache.json",
    activities_out: str | Path = "artifacts/chembl_activities.parquet",
    targets_out: str | Path = "artifacts/targets.parquet",
    edges_out: str | Path = "artifacts/evidence_edges.parquet",
    figure_out: str | Path = "figures/chembl_and_hallmark_coverage.png",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Execute Phase 2 ChEMBL join and evidence grading pipeline.

    Returns:
        Tuple of (coverage_df, activities_df, targets_df, edges_df).
    """
    comp_file = Path(compounds_parquet)
    if not comp_file.exists():
        raise FileNotFoundError(f"Compounds parquet file not found: {comp_file}")

    logger.info("Reading resolved compounds from %s...", comp_file)
    compounds_df = pd.read_parquet(comp_file)

    client = ChEMBLClient(cache_path=cache_path)

    coverage_records = []
    all_activities = []
    unique_target_ids = set()

    for _, row in compounds_df.iterrows():
        raw_name = str(row["raw_name"])
        inchikey = row.get("inchikey")
        resolved = bool(row.get("resolved", False))

        chembl_id = None
        has_chembl = False
        has_binding = False
        num_acts = 0
        num_binding = 0

        if resolved and inchikey and pd.notna(inchikey):
            logger.info("Querying ChEMBL for '%s' (InChIKey: %s)...", raw_name, inchikey)
            chembl_id = client.get_chembl_id_by_inchikey(inchikey)

            if chembl_id:
                has_chembl = True
                acts = client.get_compound_activities(chembl_id, limit=50)
                num_acts = len(acts)

                for act in acts:
                    assay_type = str(act.get("assay_type", "")).upper()
                    if assay_type == "B":
                        has_binding = True
                        num_binding += 1

                    target_id = act.get("target_chembl_id")
                    if target_id:
                        unique_target_ids.add(target_id)

                    all_activities.append(
                        {
                            "compound_name": raw_name,
                            "inchikey": inchikey,
                            "molecule_chembl_id": chembl_id,
                            **act,
                        }
                    )

        coverage_records.append(
            {
                "raw_name": raw_name,
                "inchikey": inchikey,
                "chembl_id": chembl_id,
                "has_chembl": has_chembl,
                "has_binding": has_binding,
                "num_activities": num_acts,
                "num_binding_assays": num_binding,
            }
        )

    coverage_df = pd.DataFrame(coverage_records)
    activities_df = pd.DataFrame(all_activities)

    # Fetch target details
    target_records = []
    logger.info("Fetching details for %d unique target records...", len(unique_target_ids))
    for t_id in unique_target_ids:
        t_meta = client.get_target_details(t_id)
        if t_meta:
            target_records.append(t_meta)
    targets_df = pd.DataFrame(target_records)

    # Ingest and validate curated evidence edges
    curated_file = Path(curated_edges_csv)
    if not curated_file.exists():
        raise FileNotFoundError(f"Curated evidence file not found: {curated_file}")

    logger.info("Ingesting curated evidence edges from %s...", curated_file)
    edges_raw = pd.read_csv(curated_file)

    graph = EvidenceGraph()
    validated_edges: list[dict[str, Any]] = []

    for _, e_row in edges_raw.iterrows():
        try:
            grade_enum = EvidenceGrade(str(e_row["grade"]).strip().upper())
        except ValueError as exc:
            raise ValueError(
                f"Invalid EvidenceGrade in row: {e_row.to_dict()}. Grade must be E0-E4."
            ) from exc

        hallmark_val = str(e_row.get("hallmark", "")).strip().lower()
        try:
            hallmark_enum = HallmarkSlug(hallmark_val) if hallmark_val else None
        except ValueError:
            hallmark_enum = None

        doc_ids_raw = str(e_row.get("document_ids", "")).strip()
        doc_ids = (
            [d.strip() for d in doc_ids_raw.split(";") if d.strip()]
            if doc_ids_raw and doc_ids_raw != "None"
            else []
        )

        edge = EvidenceEdge(
            compound_inchikey=str(e_row.get("inchikey", "")),
            relation=str(e_row["relation"]),
            grade=grade_enum,
            target_id=str(e_row.get("target_id", "")) if pd.notna(e_row.get("target_id")) else None,
            hallmark=hallmark_enum,
            provenance=str(e_row.get("provenance", "")),
            document_ids=doc_ids,
        )
        graph.add_edge(edge)

        validated_edges.append(
            {
                "compound_name": str(e_row["compound_name"]),
                "inchikey": edge.compound_inchikey,
                "target_id": edge.target_id,
                "target_symbol": str(e_row.get("target_symbol", "")),
                "hallmark": edge.hallmark.value if edge.hallmark else None,
                "relation": edge.relation,
                "grade": edge.grade.value,
                "provenance": edge.provenance,
                "document_ids": ";".join(edge.document_ids),
                "notes": str(e_row.get("notes", "")),
            }
        )

    edges_df = pd.DataFrame(validated_edges)

    # Save artifacts
    for p, df in [
        (activities_out, activities_df),
        (targets_out, targets_df),
        (edges_out, edges_df),
        ("artifacts/chembl_coverage.parquet", coverage_df),
    ]:
        out_p = Path(p)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(out_p, index=False)
        logger.info("Saved %d records to %s", len(df), out_p)

    # Save CSV copies
    Path("data/processed").mkdir(parents=True, exist_ok=True)
    activities_df.to_csv("data/processed/chembl_activities.csv", index=False)
    targets_df.to_csv("data/processed/targets.csv", index=False)
    edges_df.to_csv("data/processed/evidence_edges.csv", index=False)

    # Generate visualization
    fig_path = plot_chembl_and_hallmark_coverage(coverage_df, edges_df, output_path=figure_out)
    logger.info("Saved ChEMBL and hallmark coverage figure to %s", fig_path)

    # Headline metric
    total_compounds = len(coverage_df)
    has_binding_count = coverage_df["has_binding"].sum()
    has_chembl_count = coverage_df["has_chembl"].sum()
    binding_pct = (has_binding_count / total_compounds * 100) if total_compounds > 0 else 0
    chembl_pct = (has_chembl_count / total_compounds * 100) if total_compounds > 0 else 0

    logger.info("=" * 60)
    logger.info("PHASE 2 SUMMARY METRICS:")
    logger.info("Total Watchlist Compounds: %d", total_compounds)
    logger.info("Compounds in ChEMBL: %d (%.1f%%)", has_chembl_count, chembl_pct)
    logger.info("Compounds with >= 1 Binding Assay: %d (%.1f%%)", has_binding_count, binding_pct)
    logger.info("Total ChEMBL Bioactivities Retrieved: %d", len(activities_df))
    logger.info("Total Unique Targets Mapped: %d", len(targets_df))
    logger.info("Curated Evidence Edges (E0-E4): %d", len(edges_df))
    logger.info("=" * 60)

    return coverage_df, activities_df, targets_df, edges_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Phase 2 ChEMBL Join Pipeline")
    parser.add_argument("--compounds", default="artifacts/compounds.parquet")
    parser.add_argument("--curated", default="data/curated/curated_evidence.csv")
    parser.add_argument("--cache", default="data/interim/chembl_cache.json")
    args = parser.parse_args()

    run_phase2_pipeline(
        compounds_parquet=args.compounds,
        curated_edges_csv=args.curated,
        cache_path=args.cache,
    )
