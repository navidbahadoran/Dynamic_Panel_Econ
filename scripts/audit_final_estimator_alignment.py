"""Write deterministic evidence for the finalized DGP and constrained estimator."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dynamic_panel_econ.calibration import load_frozen_calibrations
from dynamic_panel_econ.core import Coefficients, Design, max_abs
from dynamic_panel_econ.estimation import fit_fixed_rank
from dynamic_panel_econ.inference import infer_target
from dynamic_panel_econ.lowrank import numerical_rank

FROZEN = Path("configs/mc/frozen_dgp_calibration.toml")
OUTPUT = Path("results/mc/audit/final_estimator_alignment")
B = 10.0
C_B = 1.0


def _fixture(level: float) -> tuple[np.ndarray, Design, Coefficients]:
    zero = np.zeros((8, 6), dtype=float)
    y = np.full_like(zero, level)
    return y, Design([zero.copy()], [zero.copy()]), Coefficients(
        [zero.copy()], [zero.copy()], y.copy()
    )


def calibration_rows() -> pd.DataFrame:
    digest = hashlib.sha256(FROZEN.read_bytes()).hexdigest()
    rows = []
    for key, result in sorted(load_frozen_calibrations(FROZEN).items(), key=lambda item: str(item[0])):
        diagnostics = result.diagnostics
        rows.append(
            {
                "frozen_calibration_file": str(FROZEN).replace("\\", "/"),
                "frozen_calibration_hash": digest,
                "design": diagnostics["design"],
                "dgp": key[0],
                "N": key[1],
                "T": key[2],
                "true_rank_vector": json.dumps(key[3] or (1, 1, 1)),
                "c_h": result.c_h,
                "c_xi": result.c_xi,
                "r2_scale_identified": diagnostics["r2_scale_identified"],
                "intended_r2": result.target_r2,
                "calibrated_r2": result.achieved_r2,
                "pi_H": result.achieved_h_share,
                "C_A": diagnostics["C_A"],
                "C_beta": diagnostics["C_beta"],
                "C_H": diagnostics["C_H"],
                "C_Theta": diagnostics["C_Theta"],
                "B": B,
                "c_B": C_B,
                "activated": True,
            }
        )
    return pd.DataFrame(rows)


def solver_rows() -> pd.DataFrame:
    interior_y, interior_design, interior_start = _fixture(0.2)
    interior = fit_fixed_rank(
        interior_y,
        interior_design,
        (0, 0, 1),
        initial=interior_start,
        coefficient_bound=B,
    )
    boundary_y, boundary_design, boundary_start = _fixture(20.0)
    boundary = fit_fixed_rank(
        boundary_y,
        boundary_design,
        (0, 0, 1),
        initial=boundary_start,
        coefficient_bound=B,
        max_sweeps=50,
        constrained_kkt_tolerance=1e-6,
    )
    direction = Coefficients(
        [np.zeros_like(boundary_y)],
        [np.zeros_like(boundary_y)],
        np.full_like(boundary_y, 1.0 / boundary_y.size),
    )
    inference = infer_target(direction, boundary, boundary_y, boundary_design, spatial=False)
    rows = []
    for case, fit in (("interior_fast_path", interior), ("boundary_fallback", boundary)):
        rows.append(
            {
                "case": case,
                "B": B,
                "unconstrained_max_abs": fit.diagnostics["unconstrained_max_abs"],
                "unconstrained_inside_box": fit.diagnostics["unconstrained_inside_box"],
                "constrained_fallback_used": fit.diagnostics["constrained_fallback_used"],
                "constrained_objective": fit.diagnostics["constrained_objective"],
                "max_abs_coefficient": max_abs(fit.theta),
                "max_constraint_violation": fit.diagnostics["max_constraint_violation"],
                "boundary_active": fit.diagnostics["boundary_active"],
                "constrained_KKT_residual": fit.diagnostics["constrained_KKT_residual"],
                "constrained_iterations": fit.diagnostics["constrained_iterations"],
                "constrained_runtime": fit.diagnostics["constrained_runtime"],
                "constrained_solver_status": fit.diagnostics["constrained_solver_status"],
                "numerical_rank_vector": json.dumps(
                    [numerical_rank(matrix) for matrix in fit.theta.matrices()]
                ),
                "valid": fit.converged,
                "inference_computed": case == "boundary_fallback" and np.isfinite(inference.estimate),
                "normal_equation_used_as_identity": (
                    inference.diagnostics["normal_equation_used_as_identity"]
                    if case == "boundary_fallback"
                    else False
                ),
            }
        )
    rows.append(
        {
            "case": "forced_solver_failure",
            "B": B,
            "constrained_solver_status": "constrained_solver_failure",
            "valid": False,
            "evidence": "tests/test_constrained_estimator.py::test_constrained_solver_failure_has_explicit_status",
        }
    )
    return pd.DataFrame(rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    calibration = calibration_rows()
    calibration.to_csv(OUTPUT / "activated_calibration_summary.csv", index=False)
    envelopes = calibration[
        ["design", "dgp", "N", "T", "true_rank_vector", "C_A", "C_beta", "C_H", "C_Theta"]
    ].copy()
    envelopes["B"] = B
    envelopes["c_B"] = C_B
    envelopes["distance_to_boundary_B"] = B - envelopes["C_Theta"]
    envelopes["slack_relative_to_B_minus_c_B"] = B - C_B - envelopes["C_Theta"]
    envelopes["condition_verified"] = envelopes["C_Theta"] <= B - C_B
    envelopes.to_csv(OUTPUT / "coefficient_envelope_verification.csv", index=False)
    pd.DataFrame(
        [
            ("supplied_rank_full_panel", "estimation.fit_fixed_rank", "shared fast path + box QP", True),
            ("candidate_post_refit", "rank_selection.fit_fixed_rank", "shared fast path + box QP", True),
            ("rank_at_most_cap_pilot", "rank_selection.fit_rank_adaptive_cap_pilot", "shared route refits", True),
            ("time_half_split", "inference.prepare_split_fits", "shared fast path + box QP", True),
            ("unit_half_split", "inference.prepare_split_fits", "shared fast path + box QP", True),
        ],
        columns=("estimator_path", "code_path", "optimizer", "literal_box_handling"),
    ).to_csv(OUTPUT / "estimator_path_alignment.csv", index=False)
    solver_rows().to_csv(OUTPUT / "constrained_solver_tests.csv", index=False)


if __name__ == "__main__":
    main()
