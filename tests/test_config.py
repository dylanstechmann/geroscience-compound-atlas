"""Tests for benchmark configuration loading and defaults."""

from pathlib import Path

from bench.config import (
    BenchConfig,
    load_bench_config,
)


def test_default_config_construction():
    """Verify BenchConfig defaults match expected values without a YAML file."""
    cfg = BenchConfig()

    assert cfg.primary_task == "senolytic_selectivity"
    assert cfg.fallback_task == "chembl_mtor_activity"
    assert cfg.split.type == "bemis_murcko"
    assert cfg.split.train_ratio == 0.80
    assert cfg.split.val_ratio == 0.10
    assert cfg.split.test_ratio == 0.10
    assert cfg.split.seeds == (42, 123, 456)
    assert cfg.features.fingerprint == "morgan"
    assert cfg.features.radius == 2
    assert cfg.features.n_bits == 2048
    assert cfg.features.include_descriptors is True
    assert cfg.model.baseline == "logistic_regression"
    assert cfg.model.contender == "hist_gradient_boosting"
    assert cfg.model.hgb_params.max_iter == 200
    assert cfg.model.hgb_params.learning_rate == 0.05
    assert cfg.model.hgb_params.min_samples_leaf == 10
    assert cfg.thresholds.pchembl_active == 6.0
    assert "auroc" in cfg.metrics
    assert "brier_score" in cfg.metrics


def test_load_bench_config_from_yaml():
    """Verify config loading from the real bench.yaml file."""
    config_path = Path(__file__).resolve().parent.parent / "configs" / "bench.yaml"
    if not config_path.exists():
        # Fall back to compiled defaults if config file missing in test env
        cfg = load_bench_config(config_path)
        assert isinstance(cfg, BenchConfig)
        return

    cfg = load_bench_config(config_path)
    assert isinstance(cfg, BenchConfig)
    assert cfg.split.type == "bemis_murcko"
    assert cfg.model.hgb_params.max_iter == 200
    assert cfg.thresholds.pchembl_active == 6.0


def test_load_bench_config_missing_file_uses_defaults():
    """Verify that a missing config file gracefully falls back to defaults."""
    cfg = load_bench_config("/nonexistent/path/bench.yaml")
    assert isinstance(cfg, BenchConfig)
    assert cfg.primary_task == "senolytic_selectivity"
    assert cfg.split.seeds == (42, 123, 456)


def test_frozen_dataclasses():
    """Verify config dataclasses are immutable (frozen=True)."""
    cfg = BenchConfig()
    try:
        cfg.primary_task = "something_else"  # type: ignore[misc]
        assert False, "Should have raised FrozenInstanceError"
    except AttributeError:
        pass  # Expected: frozen dataclass
