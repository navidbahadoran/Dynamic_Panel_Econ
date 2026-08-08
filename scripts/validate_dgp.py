"""Run deterministic calibration and DGP diagnostics, including DGP 4 truths."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from dynamic_panel_econ.calibration import calibrate_cell, pooled_r2
from dynamic_panel_econ.config import config_hash, load_config
from dynamic_panel_econ.core import fitted_values
from dynamic_panel_econ.dgp import generate_panel
from dynamic_panel_econ.lowrank import numerical_rank
from dynamic_panel_econ.monte_carlo import _params, resolve_group_gap
from dynamic_panel_econ.seeds import seed_sequence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config, gap_reports = resolve_group_gap(load_config(args.config))
    params = _params(config)
    output = Path(config["run"]["output_root"]) / "calibration" / config_hash(config)
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    truth_rows = []
    for dgp in config["run"]["dgps"]:
        for n, t in config["run"]["cells"]:
            calibration = calibrate_cell(
                int(dgp), int(n), int(t), int(config["run"]["master_seed"]),
                params=params,
                pi_h=float(config["dgp"]["pi_h"]),
                target_r2=float(config["dgp"]["target_r2"]),
                draws=int(config["dgp"]["calibration_draws"]),
            )
            row = asdict(calibration)
            row.update(row.pop("diagnostics"))
            panel = generate_panel(
                int(dgp), int(n), int(t),
                seed_sequence(config["run"]["master_seed"], "validation", dgp, n, t),
                c_h=calibration.c_h, c_xi=calibration.c_xi, params=params,
            )
            fitted = fitted_values(panel.theta0, panel.design)
            row["validation_r2"] = pooled_r2(panel.y, fitted)
            row["corr_x_current_u"] = float(np.corrcoef(panel.design.x[0].ravel(), panel.u_tilde.ravel())[0, 1])
            row["corr_x_lagged_u"] = float(np.corrcoef(panel.design.x[0].ravel(), panel.u_tilde_lag.ravel())[0, 1])
            row["true_ranks"] = json.dumps([numerical_rank(m) for m in panel.theta0.matrices()])
            rows.append(row)
            if int(dgp) == 4:
                truth_rows.append({"N": n, "T": t, **panel.diagnostics, **panel.truths})
    pd.DataFrame(rows).to_parquet(output / "dgp_validation.parquet", index=False)
    pd.DataFrame(rows).to_csv(output / "dgp_validation.csv", index=False)
    pd.DataFrame(truth_rows).to_csv(output / "dgp4_group_truths.csv", index=False)
    (output / "dgp4_group_gap_pilot.json").write_text(json.dumps(gap_reports, indent=2) + "\n", encoding="utf-8")
    print(pd.DataFrame(rows).to_string(index=False))
    if truth_rows:
        columns = [
            "N", "T", "mu_lambda_a_1", "mu_lambda_a_2", "c_a",
            "A_G1_fixed_time_true", "A_G2_fixed_time_true", "A_G2_minus_G1_fixed_time_true",
            "A_G1_time_average_true", "A_G2_time_average_true", "A_G2_minus_G1_time_average_true",
        ]
        print("\nDGP 4 exact post-scaling autoregressive group truths:\n")
        print(pd.DataFrame(truth_rows)[columns].to_string(index=False))
    print(f"\nValidation outputs: {output}")


if __name__ == "__main__":
    main()
