"""Command-line interface for Geroscience Compound Atlas."""

from __future__ import annotations

import argparse
import sys
from typing import Sequence

from atlas.lookup import main as lookup_main
from atlas.pipeline import run_phase1_pipeline
from atlas.chembl_join import run_phase2_pipeline
from viz.dashboard import build_dashboard_html


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="geroatlas",
        description="Geroscience Compound Atlas & Predictive Bench CLI",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: lookup
    lookup_parser = subparsers.add_parser(
        "lookup",
        help="Single-compound lookup (evidence card, hallmark edges, benchmark prediction)",
    )
    lookup_parser.add_argument("compound", help="Compound name, synonym, CID, InChIKey, or canonical SMILES")
    lookup_parser.add_argument("--json", action="store_true", help="Output machine-readable JSON")
    lookup_parser.add_argument("--out", type=str, help="Optional output file destination")
    lookup_parser.add_argument("--compounds", default="artifacts/compounds.parquet", help="Path to compounds table")
    lookup_parser.add_argument("--edges", default="artifacts/evidence_edges.parquet", help="Path to evidence edges")
    lookup_parser.add_argument("--coverage", default="artifacts/chembl_coverage.parquet", help="Path to ChEMBL coverage")
    lookup_parser.add_argument("--benchmark", default="artifacts/benchmark_dataset.parquet", help="Path to benchmark parquet")
    lookup_parser.add_argument("--splits", default="artifacts/splits.json", help="Path to splits json")

    # Subcommand: pipeline
    subparsers.add_parser(
        "pipeline",
        help="Run Phase 1 compound resolution and featurization",
    )

    # Subcommand: chembl-join
    subparsers.add_parser(
        "chembl-join",
        help="Run Phase 2 ChEMBL bioactivity join and evidence graph construction",
    )

    # Subcommand: dashboard
    subparsers.add_parser(
        "dashboard",
        help="Generate standalone interactive HTML dashboard",
    )

    args = parser.parse_args(argv)

    if args.command == "lookup":
        # Forward remaining args to lookup_main
        lookup_args = [args.compound]
        if args.json:
            lookup_args.append("--json")
        if args.out:
            lookup_args.extend(["--out", args.out])
        lookup_args.extend([
            "--compounds", args.compounds,
            "--edges", args.edges,
            "--coverage", args.coverage,
            "--benchmark", args.benchmark,
            "--splits", args.splits,
        ])
        return lookup_main(lookup_args)
    elif args.command == "pipeline":
        run_phase1_pipeline()
        return 0
    elif args.command == "chembl-join":
        run_phase2_pipeline()
        return 0
    elif args.command == "dashboard":
        build_dashboard_html()
        return 0
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
