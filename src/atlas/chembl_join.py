"""Phase 2 Pipeline: Join compounds to ChEMBL bioactivities and build curated evidence graph."""

import argparse
import logging
from pathlib import Path
from typing import Any

import pandas as pd

from atlas.chembl import ChEMBLClient
from atlas.evidence_review import apply_evidence_review
from atlas.graph import EvidenceGraph
from atlas.models import EvidenceEdge, EvidenceGrade, HallmarkSlug

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("atlas.chembl_join")

# Curated ChEMBL target identifiers are paired with their approved gene symbols.
# This crosswalk is intentionally explicit: an unknown ID or contradictory pair
# must be reviewed instead of silently becoming a plausible-looking edge.
CHEMBL_TARGET_SYMBOLS = {
    "CHEMBL2842": "MTOR",
    "CHEMBL2814": "SRC",
    "CHEMBL219": "EPHA2",
    "CHEMBL4203": "SIRT1",
    "CHEMBL376": "PRKAA1",
    "CHEMBL1201509": "EP300",
    "CHEMBL246": "ESR1",
    "CHEMBL444141": "BCL2",
    "CHEMBL2047": "BCL2L1",
    "CHEMBL1613": "MGAM",
    "CHEMBL429815": "NAMPT",
    "CHEMBL1926673": "NMRK1",
    "CHEMBL4126937": "PINK1",
    "CHEMBL2103833": "MAP2K1",
    "CHEMBL2029969": "SLC5A2",
}

REQUIRED_CURATED_COLUMNS = {
    "compound_name", "compound_id", "inchikey", "target_id", "target_symbol",
    "hallmark", "relation", "grade", "provenance", "document_ids", "notes",
    "source_url", "source_review_status", "evidence_basis", "intervention_components",
    "attribution_scope",
}


def _normalized_optional(value: Any) -> str | None:
    if pd.isna(value):
        return None
    normalized = str(value).strip()
    return normalized or None


def _compound_id(cid: Any) -> str | None:
    if pd.isna(cid):
        return None
    return f"pubchem:{int(float(cid))}"


def validate_curated_evidence(
    curated_df: pd.DataFrame, compounds_df: pd.DataFrame
) -> None:
    """Fail closed on dangling identities, target conflicts, and unsupported E4 rows."""
    missing_columns = REQUIRED_CURATED_COLUMNS - set(curated_df.columns)
    if missing_columns:
        raise ValueError(f"Curated evidence is missing required columns: {sorted(missing_columns)}")
    required_compound_columns = {"raw_name", "cid", "inchikey", "resolved"}
    missing_compound_columns = required_compound_columns - set(compounds_df.columns)
    if missing_compound_columns:
        raise ValueError(f"Compound snapshot is missing required columns: {sorted(missing_compound_columns)}")

    resolved = compounds_df[compounds_df["resolved"].fillna(False)].copy()
    resolved["_compound_id"] = resolved["cid"].map(_compound_id)
    if resolved["_compound_id"].isna().any() or resolved["_compound_id"].duplicated().any():
        raise ValueError("Resolved compound snapshot must have unique, nonempty PubChem CIDs")
    by_id = resolved.set_index("_compound_id")
    by_name = resolved.groupby("raw_name", dropna=False)

    for idx, row in curated_df.iterrows():
        compound_name = _normalized_optional(row["compound_name"])
        compound_id = _normalized_optional(row["compound_id"])
        inchikey = _normalized_optional(row["inchikey"])
        if not compound_name or not compound_id or not inchikey:
            raise ValueError(f"Curated row {idx} has an empty compound identity")
        if compound_id not in by_id.index:
            raise ValueError(f"Curated row {idx} has unresolved compound_id {compound_id}")
        if compound_name not in by_name.groups or len(by_name.get_group(compound_name)) != 1:
            raise ValueError(f"Curated row {idx} compound name is not an exact unique snapshot name")
        compound = by_id.loc[compound_id]
        if str(compound["raw_name"]) != compound_name:
            raise ValueError(f"Curated row {idx} compound_id does not match exact name {compound_name}")
        if str(compound["inchikey"]) != inchikey:
            raise ValueError(f"Curated row {idx} InChIKey disagrees with {compound_id}")

        target_id = _normalized_optional(row["target_id"])
        target_symbol = _normalized_optional(row["target_symbol"])
        if target_id:
            expected_symbol = CHEMBL_TARGET_SYMBOLS.get(target_id)
            if expected_symbol is None:
                raise ValueError(f"Curated row {idx} has unreviewed target ID {target_id}")
            if not target_symbol or target_symbol.upper() != expected_symbol:
                raise ValueError(
                    f"Curated row {idx} target conflict: {target_id} maps to {expected_symbol}, "
                    f"not {target_symbol!r}"
                )

        if not _normalized_optional(row["source_review_status"]):
            raise ValueError(f"Curated row {idx} has no source review status")
        grade = str(row["grade"]).strip().upper()
        source_review_status = _normalized_optional(row["source_review_status"])
        allowed_review_statuses = {
            "claim_support_unreviewed",
            "source_reviewed_with_claim_limits",
            "unverified_vendor_or_gray_market_claim",
        }
        if source_review_status not in allowed_review_statuses:
            raise ValueError(f"Curated row {idx} has an unsupported source review status")
        curator_grade = (
            _normalized_optional(row.get("curator_grade")) or grade
        ).upper()
        try:
            EvidenceGrade(grade)
            EvidenceGrade(curator_grade)
        except ValueError as exc:
            raise ValueError(f"Curated row {idx} has an invalid displayed or curator grade") from exc
        evidence_basis = _normalized_optional(row["evidence_basis"])
        if (grade == "E4" or curator_grade == "E4") and evidence_basis not in {
            "mammalian_lifespan", "reported_human_intervention_outcome"
        }:
            raise ValueError(f"Curated row {idx} uses E4 without a reported lifespan/outcome")
        expected_grade = (
            curator_grade
            if source_review_status == "source_reviewed_with_claim_limits"
            else EvidenceGrade.E0.value
        )
        if grade != expected_grade:
            raise ValueError(
                f"Curated row {idx} displayed grade {grade} conflicts with review status; "
                f"expected {expected_grade} from curator grade {curator_grade}"
            )

        source_url = _normalized_optional(row["source_url"])
        if not source_url and evidence_basis != "unverified_vendor_claim":
            raise ValueError(f"Curated row {idx} has no source URL and is not explicitly marked as an unverified claim")

        document_ids = _normalized_optional(row["document_ids"]) or ""
        source_url = source_url or ""
        for doc_id in (part.strip() for part in document_ids.split(";")):
            if doc_id.startswith("PMID:"):
                expected_url = f"https://pubmed.ncbi.nlm.nih.gov/{doc_id.removeprefix('PMID:')}"
            elif doc_id.startswith("CHEMBL"):
                expected_url = f"https://www.ebi.ac.uk/chembl/explore/assay/{doc_id}"
            else:
                continue
            if expected_url not in source_url:
                raise ValueError(f"Curated row {idx} is missing a direct source URL for {doc_id}")


def build_curated_edges_frame(edges_raw: pd.DataFrame) -> pd.DataFrame:
    """Validate and materialize local curated edges without contacting ChEMBL."""
    graph = EvidenceGraph()
    validated_edges: list[dict[str, Any]] = []
    for _, e_row in edges_raw.iterrows():
        source_review_status = _normalized_optional(e_row.get("source_review_status")) or ""
        try:
            curator_grade = EvidenceGrade(
                (_normalized_optional(e_row.get("curator_grade")) or str(e_row["grade"])).strip().upper()
            )
        except ValueError as exc:
            raise ValueError(
                f"Invalid EvidenceGrade in row: {e_row.to_dict()}. Grade must be E0-E4."
            ) from exc

        # Curator-assigned grades are not promoted until the source-to-claim
        # relationship has been reviewed. Keep the asserted grade separately
        # so the record remains auditable without presenting it as verified.
        display_grade = (
            curator_grade
            if source_review_status == "source_reviewed_with_claim_limits"
            else EvidenceGrade.E0
        )

        hallmark_val = (_normalized_optional(e_row.get("hallmark")) or "").lower()
        try:
            hallmark_enum = HallmarkSlug(hallmark_val) if hallmark_val else None
        except ValueError:
            hallmark_enum = None

        doc_ids_raw = _normalized_optional(e_row.get("document_ids")) or ""
        doc_ids = (
            [d.strip() for d in doc_ids_raw.split(";") if d.strip()]
            if doc_ids_raw and doc_ids_raw != "None"
            else []
        )
        edge = EvidenceEdge(
            compound_inchikey=str(e_row.get("inchikey", "")),
            relation=str(e_row["relation"]),
            grade=display_grade,
            target_id=_normalized_optional(e_row.get("target_id")),
            hallmark=hallmark_enum,
            provenance=str(e_row.get("provenance", "")),
            document_ids=doc_ids,
        )
        graph.add_edge(edge)
        validated_edges.append(
            {
                "compound_name": str(e_row["compound_name"]),
                "compound_id": str(e_row["compound_id"]),
                "inchikey": edge.compound_inchikey,
                "target_id": edge.target_id,
                "target_symbol": _normalized_optional(e_row.get("target_symbol")),
                "hallmark": edge.hallmark.value if edge.hallmark else None,
                "relation": edge.relation,
                "grade": edge.grade.value,
                "curator_grade": curator_grade.value,
                "provenance": edge.provenance,
                "document_ids": ";".join(edge.document_ids),
                "notes": _normalized_optional(e_row.get("notes")) or "",
                "source_url": _normalized_optional(e_row.get("source_url")) or "",
                "source_title": _normalized_optional(e_row.get("source_title")),
                "species": _normalized_optional(e_row.get("species")),
                "sex": _normalized_optional(e_row.get("sex")),
                "study_design": _normalized_optional(e_row.get("study_design")),
                "endpoint": _normalized_optional(e_row.get("endpoint")),
                "comparator": _normalized_optional(e_row.get("comparator")),
                "effect_estimate": _normalized_optional(e_row.get("effect_estimate")),
                "uncertainty": _normalized_optional(e_row.get("uncertainty")),
                "source_locator": _normalized_optional(e_row.get("source_locator")),
                "evidence_basis": _normalized_optional(e_row.get("evidence_basis")) or "",
                "source_review_status": source_review_status,
                "intervention_components": _normalized_optional(e_row.get("intervention_components")) or "",
                "attribution_scope": _normalized_optional(e_row.get("attribution_scope")) or "",
            }
        )
    return apply_evidence_review(pd.DataFrame(validated_edges))


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

    curated_file = Path(curated_edges_csv)
    if not curated_file.exists():
        raise FileNotFoundError(f"Curated evidence file not found: {curated_file}")
    logger.info("Validating curated evidence identities and source metadata...")
    edges_raw = pd.read_csv(curated_file)
    validate_curated_evidence(edges_raw, compounds_df)

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
    logger.info("Ingesting curated evidence edges from %s...", curated_file)

    edges_df = build_curated_edges_frame(edges_raw)

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
    parser.add_argument("--figure", default="figures/chembl_and_hallmark_coverage.png")
    parser.add_argument(
        "--curated-only", action="store_true",
        help="rebuild only evidence_edges.parquet from local sources; never initializes ChEMBL",
    )
    parser.add_argument("--edges-out", default="artifacts/evidence_edges.parquet")
    args = parser.parse_args()

    if args.curated_only:
        compounds_df = pd.read_parquet(args.compounds)
        curated_df = pd.read_csv(args.curated)
        validate_curated_evidence(curated_df, compounds_df)
        edges_df = build_curated_edges_frame(curated_df)
        edges_path = Path(args.edges_out)
        edges_path.parent.mkdir(parents=True, exist_ok=True)
        edges_df.to_parquet(edges_path, index=False)
        logger.info("Saved %d local curated edges to %s", len(edges_df), edges_path)
        raise SystemExit(0)

    cov_df, acts_df, tgts_df, edg_df = run_phase2_pipeline(
        compounds_parquet=args.compounds,
        curated_edges_csv=args.curated,
        cache_path=args.cache,
    )

    # Generate visualization (kept out of the data pipeline to avoid viz dependency)
    from viz.hallmarks import plot_chembl_and_hallmark_coverage

    fig_path = plot_chembl_and_hallmark_coverage(cov_df, edg_df, output_path=args.figure)
    logger.info("Saved ChEMBL and hallmark coverage figure to %s", fig_path)
