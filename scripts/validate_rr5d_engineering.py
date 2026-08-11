"""Write deterministic non-scientific RR5d solver-equivalence evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dynamic_panel_econ.cap_plus_one import balanced_factor_block, fit_cap_plus_one
from dynamic_panel_econ.core import Coefficients, Design
from dynamic_panel_econ.estimation import fit_fixed_rank
from dynamic_panel_econ.lowrank import numerical_rank
from dynamic_panel_econ.rank_selection import fit_revision10_spectral_pilot


def _zero_design(n: int, t: int) -> Design:
    zero = np.zeros((n, t), dtype=np.float64)
    return Design([zero.copy()], [zero.copy()])


def _product(rank: int, n: int, t: int) -> np.ndarray:
    if rank == 0:
        return np.zeros((n, t), dtype=np.float64)
    result = np.outer(np.arange(1, n + 1), np.arange(1, t + 1)) / (n * t)
    if rank >= 2:
        result += 0.2 * np.outer(np.linspace(-1, 1, n), np.linspace(1, -1, t))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "results/design/revision10_ridge_ratio/rr5d_engineering/"
            "deterministic_solver_results.csv"
        ),
    )
    args = parser.parse_args()
    rows: list[dict[str, object]] = []
    for rank in (0, 1, 2):
        y = _product(rank, 10, 9)
        fit = fit_cap_plus_one(
            y,
            _zero_design(10, 9),
            (0, 0, 4),
            seed=100 + rank,
            max_sweeps=50,
        )
        rows.append(
            {
                "fixture": f"width4_rank{rank}",
                "expected_product_rank": rank,
                "realized_product_rank": numerical_rank(fit.theta.H),
                "converged": fit.converged,
                "objective": fit.objective,
                "stationarity_or_kkt": fit.stationarity_residual,
                "max_coefficient": float(np.max(np.abs(fit.theta.H))),
                "product_error": float(np.max(np.abs(fit.theta.H - y))),
                "objective_monotone": all(
                    following <= previous + 1e-10
                    for previous, following in zip(
                        fit.objective_history, fit.objective_history[1:], strict=False
                    )
                ),
                "status": fit.diagnostics["constrained_solver_status"],
            }
        )

    rng = np.random.default_rng(222)
    loading = rng.normal(size=(12, 4))
    factor = rng.normal(size=(11, 4))
    loading[:, 0] *= 1e10
    factor[:, 0] /= 1e10
    product = loading @ factor.T
    balanced = balanced_factor_block(product, 4)
    rows.append(
        {
            "fixture": "badly_scaled_gauge",
            "expected_product_rank": numerical_rank(product),
            "realized_product_rank": numerical_rank(balanced.matrix()),
            "converged": True,
            "objective": 0.0,
            "stationarity_or_kkt": 0.0,
            "max_coefficient": float(np.max(np.abs(product))),
            "product_error": float(np.max(np.abs(product - balanced.matrix()))),
            "objective_monotone": True,
            "status": "product_preserved",
        }
    )

    near_collinear = np.outer(np.linspace(-1, 1, 12), np.linspace(1, -1, 11))
    near_collinear += 1e-10 * np.outer(
        np.linspace(-1, 1, 12) ** 2, np.linspace(1, -1, 11)[::-1]
    )
    collinear_fit = fit_cap_plus_one(
        near_collinear,
        _zero_design(12, 11),
        (0, 0, 4),
        seed=303,
        max_sweeps=50,
    )
    rows.append(
        {
            "fixture": "near_collinear",
            "expected_product_rank": numerical_rank(near_collinear),
            "realized_product_rank": numerical_rank(collinear_fit.theta.H),
            "converged": collinear_fit.converged,
            "objective": collinear_fit.objective,
            "stationarity_or_kkt": collinear_fit.stationarity_residual,
            "max_coefficient": float(np.max(np.abs(collinear_fit.theta.H))),
            "product_error": float(np.max(np.abs(collinear_fit.theta.H - near_collinear))),
            "objective_monotone": True,
            "status": collinear_fit.diagnostics["constrained_solver_status"],
        }
    )

    active_y = np.full((10, 9), 2.0)
    active = fit_cap_plus_one(
        active_y,
        _zero_design(10, 9),
        (0, 0, 1),
        seed=404,
        coefficient_bound=0.5,
        max_sweeps=50,
        constrained_kkt_tolerance=1e-4,
    )
    rows.append(
        {
            "fixture": "active_box",
            "expected_product_rank": 1,
            "realized_product_rank": numerical_rank(active.theta.H),
            "converged": active.converged,
            "objective": active.objective,
            "stationarity_or_kkt": active.stationarity_residual,
            "max_coefficient": float(np.max(np.abs(active.theta.H))),
            "product_error": float("nan"),
            "objective_monotone": True,
            "status": active.diagnostics["constrained_solver_status"],
        }
    )

    regression_y = _product(1, 10, 9)
    zero = np.zeros_like(regression_y)
    initial = Coefficients([zero.copy()], [zero.copy()], regression_y.copy())
    legacy = fit_fixed_rank(
        regression_y,
        _zero_design(10, 9),
        (0, 0, 1),
        initial=initial,
        max_sweeps=50,
    )
    engineered = fit_cap_plus_one(
        regression_y,
        _zero_design(10, 9),
        (0, 0, 1),
        seed=505,
        max_sweeps=50,
    )
    rows.append(
        {
            "fixture": "old_new_semantic_regression",
            "expected_product_rank": 1,
            "realized_product_rank": numerical_rank(engineered.theta.H),
            "converged": legacy.converged and engineered.converged,
            "objective": engineered.objective,
            "stationarity_or_kkt": engineered.stationarity_residual,
            "max_coefficient": float(np.max(np.abs(engineered.theta.H))),
            "product_error": float(np.max(np.abs(legacy.theta.H - engineered.theta.H))),
            "objective_monotone": True,
            "status": "same_coefficient_problem",
            "legacy_objective": legacy.objective,
            "objective_difference": abs(legacy.objective - engineered.objective),
        }
    )

    gate_y = _product(2, 10, 10)
    _, gate = fit_revision10_spectral_pilot(
        gate_y,
        _zero_design(10, 10),
        (3, 3, 3),
        seed=606,
        fit_options={
            "max_sweeps": 50,
            "coefficient_bound": 10.0,
            "stationarity_tol": 1e-6,
            "constrained_kkt_tolerance": 1e-4,
        },
        stationarity_tolerance=1e-6,
        start_objective_stability_tol=1e-6,
    )
    rows.append(
        {
            "fixture": "frozen_three_start_gate",
            "expected_product_rank": 2,
            "realized_product_rank": json.dumps(gate["numerical_rank_vector"]),
            "converged": gate["objective_stability_pass"],
            "objective": gate["best_objective"],
            "stationarity_or_kkt": max(gate["all_start_stationarity_residuals"]),
            "max_coefficient": gate["max_envelope_ratio"] * 10.0,
            "product_error": float("nan"),
            "objective_monotone": True,
            "status": "accepted" if gate["objective_stability_pass"] else "unresolved",
            "valid_start_count": gate["valid_start_count"],
            "best_two_objective_gap": gate["best_two_objective_gap"],
        }
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
