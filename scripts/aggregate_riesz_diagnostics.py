"""Aggregate the pre-production Riesz diagnostic run."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from dynamic_panel_econ.config import config_hash, load_config
from dynamic_panel_econ.diagnostics import aggregate_riesz_diagnostic, write_riesz_diagnostic
from dynamic_panel_econ.monte_carlo import resolve_group_gap


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/mc/riesz_diagnostic.toml")
    parser.add_argument("--run-root")
    parser.add_argument("--output", default="results/mc/diagnostics")
    args = parser.parse_args()
    config, _ = resolve_group_gap(load_config(args.config))
    root = (
        Path(args.run_root)
        if args.run_root
        else Path(config["run"]["output_root"]) / config["run"]["name"] / config_hash(config)
    )
    files = sorted((root / "raw").glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no raw chunks under {root / 'raw'}")
    frames = [pd.read_parquet(path) for path in files]
    raw = pd.concat([frame for frame in frames if len(frame.columns)], ignore_index=True)
    replications, events, summary = aggregate_riesz_diagnostic(
        raw,
        dgps=config["run"]["dgps"],
        cells=config["run"]["cells"],
        targets=config["inference"]["targets"],
        replications=int(config["run"]["replications"]),
        target_rayleigh_floor=float(config["inference"]["riesz_target_rayleigh_floor"]),
    )
    paths = write_riesz_diagnostic(replications, events, summary, Path(args.output))
    print(summary.to_string(index=False))
    print("\n".join(str(path) for path in paths))


if __name__ == "__main__":
    main()
