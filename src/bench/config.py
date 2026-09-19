"""Load benchmark configuration from YAML config files.

Provides a single entry point for all configurable parameters,
eliminating hardcoded magic numbers from the training and data pipelines.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)

_CONFIG_DIR = Path(__file__).resolve().parent.parent.parent / "configs"


@dataclass(frozen=True)
class SplitConfig:
    type: str = "bemis_murcko"
    train_ratio: float = 0.80
    val_ratio: float = 0.10
    test_ratio: float = 0.10
    seeds: tuple[int, ...] = (42, 123, 456)


@dataclass(frozen=True)
class FeatureConfig:
    fingerprint: str = "morgan"
    radius: int = 2
    n_bits: int = 2048
    include_descriptors: bool = True


@dataclass(frozen=True)
class HGBParams:
    max_iter: int = 200
    learning_rate: float = 0.05
    min_samples_leaf: int = 10


@dataclass(frozen=True)
class ModelConfig:
    baseline: str = "logistic_regression"
    contender: str = "hist_gradient_boosting"
    hgb_params: HGBParams = field(default_factory=HGBParams)


@dataclass(frozen=True)
class ThresholdConfig:
    pchembl_active: float = 6.0
    senolytic_selectivity_ratio: float = 2.0


@dataclass(frozen=True)
class BenchConfig:
    """Complete benchmark configuration parsed from configs/bench.yaml."""

    primary_task: str = "senolytic_selectivity"
    fallback_task: str = "chembl_mtor_activity"
    split: SplitConfig = field(default_factory=SplitConfig)
    features: FeatureConfig = field(default_factory=FeatureConfig)
    model: ModelConfig = field(default_factory=ModelConfig)
    thresholds: ThresholdConfig = field(default_factory=ThresholdConfig)
    metrics: tuple[str, ...] = ("auroc", "auprc", "recall_at_fpr_0.05", "brier_score")


def load_bench_config(config_path: str | Path | None = None) -> BenchConfig:
    """Load benchmark configuration from YAML file.

    Falls back to compiled defaults if config file is missing.
    """
    if config_path is None:
        config_path = _CONFIG_DIR / "bench.yaml"
    else:
        config_path = Path(config_path)

    if not config_path.exists():
        logger.warning("Config file not found at %s, using compiled defaults.", config_path)
        return BenchConfig()

    with open(config_path, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    bench_raw = raw.get("benchmark", {})

    split_raw = bench_raw.get("split", {})
    split_cfg = SplitConfig(
        type=split_raw.get("type", "bemis_murcko"),
        train_ratio=split_raw.get("train_ratio", 0.80),
        val_ratio=split_raw.get("val_ratio", 0.10),
        test_ratio=split_raw.get("test_ratio", 0.10),
        seeds=tuple(split_raw.get("seeds", [42, 123, 456])),
    )

    feat_raw = bench_raw.get("features", {})
    feat_cfg = FeatureConfig(
        fingerprint=feat_raw.get("fingerprint", "morgan"),
        radius=feat_raw.get("radius", 2),
        n_bits=feat_raw.get("n_bits", 2048),
        include_descriptors=feat_raw.get("include_descriptors", True),
    )

    hgb_raw = bench_raw.get("model", {}).get("hgb_params", {})
    hgb_cfg = HGBParams(
        max_iter=hgb_raw.get("max_iter", 200),
        learning_rate=hgb_raw.get("learning_rate", 0.05),
        min_samples_leaf=hgb_raw.get("min_samples_leaf", 10),
    )

    model_raw = bench_raw.get("model", {})
    model_cfg = ModelConfig(
        baseline=model_raw.get("baseline", "logistic_regression"),
        contender=model_raw.get("contender", "hist_gradient_boosting"),
        hgb_params=hgb_cfg,
    )

    threshold_raw = bench_raw.get("thresholds", {})
    threshold_cfg = ThresholdConfig(
        pchembl_active=threshold_raw.get("pchembl_active", 6.0),
        senolytic_selectivity_ratio=threshold_raw.get("senolytic_selectivity_ratio", 2.0),
    )

    metrics_raw = bench_raw.get("metrics", ["auroc", "auprc", "recall_at_fpr_0.05", "brier_score"])

    config = BenchConfig(
        primary_task=bench_raw.get("primary_task", "senolytic_selectivity"),
        fallback_task=bench_raw.get("fallback_task", "chembl_mtor_activity"),
        split=split_cfg,
        features=feat_cfg,
        model=model_cfg,
        thresholds=threshold_cfg,
        metrics=tuple(metrics_raw),
    )

    logger.info("Loaded benchmark config from %s", config_path)
    return config
