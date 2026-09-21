"""Entry point for `python -m gen`."""

import logging

from gen.pipeline import run_generator_pipeline

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

if __name__ == "__main__":
    run_generator_pipeline()
