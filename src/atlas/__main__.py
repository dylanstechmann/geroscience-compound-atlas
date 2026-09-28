"""Entry point for `python -m atlas`."""

import sys
from atlas.pipeline import run_phase1_pipeline

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] in ("lookup", "--help", "-h"):
        from atlas.cli import main
        main()
    else:
        run_phase1_pipeline()
