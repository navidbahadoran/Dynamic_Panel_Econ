"""Aggregate the frozen medium baseline and rank-stress diagnostics."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from dynamic_panel_econ.config import config_hash, load_config
from dynamic_panel_econ.medium_diagnostics import (
    aggregate_rank_medium,
    aggregate_rank_sensitivity_medium,
    aggregate_riesz_medium,
    aggregate_target_support_medium,
    attach_rank_stress_calibration,
    complete_rank_replications,
    complete_target_replications,
    write_medium_tables,
)
from dynamic_panel_econ.monte_carlo import resolve_group_gap


def _root(config_path: str) -> tuple[dict, Path]:
    config, _ = resolve_group_gap(load_config(config_path))
    root = Path(config["run"]["output_root"]) / config["run"]["name"] / config_hash(config)
    return config, root


def _tree(root: Path, subdir: str) -> pd.DataFrame:
    files = sorted((root / subdir).glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no {subdir} chunks under {root}")
    frames = [pd.read_parquet(path) for path in files]
    usable = [frame for frame in frames if len(frame.columns)]
    return pd.concat(usable, ignore_index=True) if usable else pd.DataFrame()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-config", default="configs/mc/medium_preproduction.toml")
    parser.add_argument("--stress-config", default="configs/mc/rank_stress_medium.toml")
    parser.add_argument("--output", default="results/mc/diagnostics/medium")
    args = parser.parse_args()

    baseline_config, baseline_root = _root(args.baseline_config)
    baseline_rank = _tree(baseline_root, "rank")
    baseline_raw = _tree(baseline_root, "raw")
    baseline_replications = complete_rank_replications(
        baseline_rank, baseline_raw, baseline_config, ["[1, 1, 1]"]
    )
    rank = aggregate_rank_medium(baseline_replications)
    sensitivity = aggregate_rank_sensitivity_medium(baseline_replications)
    target_replications = complete_target_replications(baseline_raw, baseline_config)
    support = aggregate_target_support_medium(target_replications)
    riesz = aggregate_riesz_medium(target_replications)

    stress_config, stress_root = _root(args.stress_config)
    stress_rank = _tree(stress_root, "rank")
    stress_raw = _tree(stress_root, "raw")
    true_ranks = [str(list(map(int, ranks))) for ranks in stress_config["rank_stress"]["true_rank_vectors"]]
    stress_replications = complete_rank_replications(
        stress_rank, stress_raw, stress_config, true_ranks
    )
    stress = aggregate_rank_medium(stress_replications).rename(
        columns={"selected_all_zero_rate": "selected_all_zero_rate"}
    )
    stress = attach_rank_stress_calibration(
        stress, pd.read_parquet(stress_root / "calibration.parquet")
    )

    paths = write_medium_tables(
        rank, support, riesz, sensitivity, stress, Path(args.output)
    )
    print(f"baseline_root={baseline_root}")
    print(f"stress_root={stress_root}")
    print("\n".join(str(path) for path in paths))


if __name__ == "__main__":
    main()
