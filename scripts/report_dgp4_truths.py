"""Summarize exact post-stability DGP 4 group truths in a deterministic pilot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dynamic_panel_econ.config import config_hash, load_config
from dynamic_panel_econ.dgp import DGP4_TRUTH_NAMES, generate_panel
from dynamic_panel_econ.monte_carlo import _params
from dynamic_panel_econ.seeds import seed_sequence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--draws", type=int, default=100)
    args = parser.parse_args()
    config = load_config(args.config)
    params = _params(config)
    rows = []
    for n, t in config["run"]["cells"]:
        draws = []
        for replication in range(args.draws):
            panel = generate_panel(
                4,
                int(n),
                int(t),
                seed_sequence(
                    config["run"]["master_seed"],
                    "dgp4_group_report",
                    n,
                    t,
                    replication,
                ),
                params=params,
            )
            draws.append({**panel.truths, "c_a": panel.diagnostics["c_a"]})
        record = {"N": int(n), "T": int(t), "draws": args.draws}
        for name in [*DGP4_TRUTH_NAMES, "c_a"]:
            values = np.array([draw[name] for draw in draws])
            for statistic, value in (
                ("mean", values.mean()),
                ("sd", values.std(ddof=1)),
                ("min", values.min()),
                ("max", values.max()),
            ):
                record[f"{name}_{statistic}"] = float(value)
        rows.append(record)
    destination = (
        Path(config["run"]["output_root"])
        / "calibration"
        / f"dgp4_truth_pilot_{config_hash(config)}"
    )
    destination.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(destination / "dgp4_truth_summary.csv", index=False)
    (destination / "dgp4_truth_summary.json").write_text(
        json.dumps(rows, indent=2) + "\n", encoding="utf-8"
    )
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"\nDGP 4 truth pilot: {destination}")


if __name__ == "__main__":
    main()
