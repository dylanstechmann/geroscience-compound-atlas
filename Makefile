.PHONY: setup data train report serve test lint format clean

PYTHON ?= python

setup:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -e ".[dev]"

data:
	$(PYTHON) -m atlas.pipeline
	$(PYTHON) -m atlas.chembl_join

train:
	$(PYTHON) -m bench.train

report:
	$(PYTHON) -m viz.dashboard

serve:
	$(PYTHON) -m http.server --directory artifacts 8000

test:
	pytest tests/ -v

lint:
	ruff check .

format:
	ruff format .

clean:
	rm -rf build/ dist/ *.egg-info .pytest_cache .ruff_cache
