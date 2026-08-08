"""Aggregate the diagnostic medium rank pilot."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from dynamic_panel_econ.config import config_hash, load_config
from dynamic_panel_econ.diagnostics import aggregate_rank_pilot, write_rank_pilot
from dynamic_panel_econ.monte_carlo import resolve_group_gap


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/mc/rank_pilot_diagnostic.toml")
    parser.add_argument("--run-root")
    parser.add_argument("--output", default="results/mc/diagnostics")
    args = parser.parse_args()
    config, _ = resolve_group_gap(load_config(args.config))
    root = (
        Path(args.run_root)
        if args.run_root
        else Path(config["run"]["output_root"]) / config["run"]["name"] / config_hash(config)
    )
    files = sorted((root / "rank").glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no rank chunks under {root / 'rank'}")
    rank_rows = pd.concat([pd.read_parquet(path) for path in files], ignore_index=True)
    replications, summary, sensitivity = aggregate_rank_pilot(
        rank_rows,
        dgps=config["run"]["dgps"],
        cells=config["run"]["cells"],
        true_rank_vectors=config["rank_stress"]["true_rank_vectors"],
        replications=int(config["run"]["replications"]),
    )
    paths = write_rank_pilot(replications, summary, sensitivity, Path(args.output))
    print(summary.to_string(index=False))
    print("\n".join(str(path) for path in paths))


if __name__ == "__main__":
    main()
