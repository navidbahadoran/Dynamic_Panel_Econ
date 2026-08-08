"""Deterministic, production-separated DGP calibration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.optimize import brentq

from .dgp import DGPParameters, _draw_rank_stress_raw, _draw_raw, _RawDraw, coefficient_envelopes
from .seeds import rng_for


@dataclass(frozen=True, slots=True)
class CalibrationResult:
    dgp: int
    n: int
    t: int
    c_h: float
    c_xi: float
    target_h_share: float
    achieved_h_share: float
    target_r2: float | None
    achieved_r2: float
    diagnostics: dict[str, Any]


class CalibrationFeasibilityError(RuntimeError):
    """Raised when the requested pooled-R2 target has no positive calibration root."""

    def __init__(self, dgp: int, n: int, t: int, target: float, lower_bound: float) -> None:
        self.dgp, self.n, self.t = dgp, n, t
        self.target, self.lower_bound = target, lower_bound
        super().__init__(
            f"DGP {dgp}, N={n}, T={t}: target pooled R2={target:.6f} is infeasible; "
            f"observed positive-c_xi lower bound is approximately {lower_bound:.6f}"
        )


def pooled_r2(y: np.ndarray, fitted: np.ndarray) -> float:
    denominator = float(np.sum((y - y.mean()) ** 2))
    if denominator <= 0.0:
        return float("nan")
    return 1.0 - float(np.sum((y - fitted) ** 2)) / denominator


def _calibrate_raws(
    dgp: int,
    n: int,
    t: int,
    raws: list[_RawDraw],
    *,
    params: DGPParameters,
    pi_h: float = 0.30,
    target_r2: float = 0.65,
    true_rank_vector: tuple[int, int, int] = (1, 1, 1),
    allow_infeasible_diagnostic: bool = False,
) -> CalibrationResult:
    """Calibrate a fixed collection of baseline or actual stress-design raw draws."""

    draws = len(raws)
    start = params.burn_in
    observed = slice(start, start + t)
    var_u = float(np.mean([np.var(raw.u_tilde[:, observed]) for raw in raws]))
    var_h = float(np.mean([np.var(raw.h_raw[:, observed]) for raw in raws]))
    if var_u <= 0.0 or var_h <= 0.0:
        raise RuntimeError("nonpositive calibration variance")
    c_h = float(np.sqrt((pi_h / (1.0 - pi_h)) * var_u / var_h))

    # Conditional on each raw draw, y(c_xi)=y_base+c_xi*y_scale.  Compute both
    # recursions once so bracketing does not rerun the dynamic simulation.
    affine_components: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for raw in raws:
        length = raw.x.shape[1]
        y_base = np.empty((n, length), dtype=np.float64)
        y_scale = np.empty((n, length), dtype=np.float64)
        previous_base = np.zeros(n, dtype=np.float64)
        previous_scale = np.zeros(n, dtype=np.float64)
        shock_scale = c_h * raw.h_raw + raw.u_tilde
        for column in range(length):
            current_base = raw.a[:, column] * previous_base + raw.beta[:, column] * raw.x[:, column]
            current_scale = raw.a[:, column] * previous_scale + shock_scale[:, column]
            y_base[:, column] = current_base
            y_scale[:, column] = current_scale
            previous_base, previous_scale = current_base, current_scale
        affine_components.append(
            (y_base[:, observed], y_scale[:, observed], raw.u_tilde[:, observed])
        )

    def average_r2(c_xi: float) -> float:
        values = []
        for y_base, y_scale, primitive_u in affine_components:
            y = y_base + c_xi * y_scale
            fitted = y - c_xi * primitive_u
            values.append(pooled_r2(y, fitted))
        return float(np.mean(values))

    r2_scale_identified = true_rank_vector[1] > 0
    bracket = None
    if r2_scale_identified:
        def root(c_xi: float) -> float:
            return average_r2(c_xi) - target_r2

        grid = np.geomspace(1e-4, 1e3, 80)
        r2_grid = [average_r2(float(value)) for value in grid]
        evaluations = [value - target_r2 for value in r2_grid]
        minimum_grid_r2 = float(min(r2_grid))
        large_c_xi_floor = float(average_r2(1e8))
        for left, right, f_left, f_right in zip(
            grid[:-1], grid[1:], evaluations[:-1], evaluations[1:], strict=True
        ):
            if f_left == 0.0 or f_left * f_right < 0.0:
                bracket = (float(left), float(right))
                break
        target_feasible: bool | None = bracket is not None
        if bracket is None and not allow_infeasible_diagnostic:
            raise CalibrationFeasibilityError(
                dgp, n, t, target_r2, min(minimum_grid_r2, large_c_xi_floor)
            )
        c_xi = (
            float(brentq(root, *bracket, xtol=1e-10, rtol=1e-10))
            if bracket is not None
            else 1.0
        )
        requested_r2: float | None = target_r2
    else:
        c_xi = 1.0
        induced = average_r2(c_xi)
        minimum_grid_r2 = induced
        large_c_xi_floor = induced
        target_feasible = None
        requested_r2 = None
    achieved_r2 = average_r2(c_xi)
    achieved_h_share = (c_h * c_h * var_h) / (c_h * c_h * var_h + var_u)
    coefficient_summary = {
        "mean_a": float(np.mean([raw.a[:, observed].mean() for raw in raws])),
        "sd_a": float(np.std(np.concatenate([raw.a[:, observed].ravel() for raw in raws]), ddof=1)),
        "mean_b": float(np.mean([raw.beta[:, observed].mean() for raw in raws])),
        "sd_b": float(np.std(np.concatenate([raw.beta[:, observed].ravel() for raw in raws]), ddof=1)),
        "min_a": float(min(raw.a[:, observed].min() for raw in raws)),
        "max_a": float(max(raw.a[:, observed].max() for raw in raws)),
        "min_b": float(min(raw.beta[:, observed].min() for raw in raws)),
        "max_b": float(max(raw.beta[:, observed].max() for raw in raws)),
        "mean_c_a": float(np.mean([raw.c_a for raw in raws])),
        "max_pre_scaling_abs_a": float(max(raw.pre_max_a for raw in raws)),
        "max_final_abs_a": float(max(np.max(np.abs(raw.a)) for raw in raws)),
        "var_u_tilde": var_u,
        "var_h_raw": var_h,
        "calibration_draws": draws,
        "minimum_grid_r2": minimum_grid_r2,
        "large_c_xi_r2_floor": large_c_xi_floor,
        "target_r2_feasible": target_feasible,
        "r2_scale_identified": r2_scale_identified,
        "requested_r2": requested_r2,
        "c_xi_normalization": 1.0 if not r2_scale_identified else None,
        "bracket_lower": bracket[0] if bracket is not None else float("nan"),
        "bracket_upper": bracket[1] if bracket is not None else float("nan"),
        "root_residual": (
            achieved_r2 - target_r2 if r2_scale_identified else None
        ),
        "infeasible_c_xi_normalization": (
            1.0 if r2_scale_identified and bracket is None else None
        ),
        "true_rank_vector": true_rank_vector,
    }
    envelopes = coefficient_envelopes(
        dgp,
        params,
        c_h=c_h,
        c_xi=c_xi,
        ranks=true_rank_vector,
    )
    coefficient_summary.update(
        {
            "theoretical_max_abs_A": envelopes["A"],
            "theoretical_max_abs_B": envelopes["B"],
            "theoretical_max_abs_H": envelopes["H"],
            "theoretical_coefficient_envelope": envelopes["all"],
            "realized_calibration_max_abs_A": float(
                max(np.max(np.abs(raw.a[:, observed])) for raw in raws)
            ),
            "realized_calibration_max_abs_B": float(
                max(np.max(np.abs(raw.beta[:, observed])) for raw in raws)
            ),
            "realized_calibration_max_abs_H": float(
                max(np.max(np.abs(c_xi * c_h * raw.h_raw[:, observed])) for raw in raws)
            ),
        }
    )
    return CalibrationResult(
        dgp=dgp,
        n=n,
        t=t,
        c_h=c_h,
        c_xi=c_xi,
        target_h_share=pi_h,
        achieved_h_share=float(achieved_h_share),
        target_r2=requested_r2,
        achieved_r2=achieved_r2,
        diagnostics=coefficient_summary,
    )


def calibrate_cell(
    dgp: int,
    n: int,
    t: int,
    master_seed: int,
    *,
    params: DGPParameters | None = None,
    pi_h: float = 0.30,
    target_r2: float = 0.65,
    draws: int = 3,
) -> CalibrationResult:
    """Calibrate ``c_H`` and ``c_xi`` with deterministic common random numbers."""

    params = params or DGPParameters()
    raws = [
        _draw_raw(dgp, n, t, rng_for(master_seed, "calibration", dgp, n, t, j), params)
        for j in range(draws)
    ]
    return _calibrate_raws(
        dgp,
        n,
        t,
        raws,
        params=params,
        pi_h=pi_h,
        target_r2=target_r2,
    )


def calibrate_rank_stress_cell(
    dgp: int,
    n: int,
    t: int,
    true_rank_vector: tuple[int, int, int],
    master_seed: int,
    *,
    component_strengths: tuple[float, ...],
    params: DGPParameters | None = None,
    pi_h: float = 0.30,
    target_r2: float = 0.65,
    draws: int = 3,
    allow_infeasible_diagnostic: bool = False,
) -> CalibrationResult:
    """Calibrate using the actual rescaled rank-stress coefficient matrices."""

    params = params or DGPParameters()
    raws = [
        _draw_rank_stress_raw(
            dgp,
            n,
            t,
            true_rank_vector,
            rng_for(
                master_seed,
                "rank_stress_calibration",
                dgp,
                n,
                t,
                true_rank_vector,
                j,
            ),
            params,
            component_strengths,
        )
        for j in range(draws)
    ]
    return _calibrate_raws(
        dgp,
        n,
        t,
        raws,
        params=params,
        pi_h=pi_h,
        target_r2=target_r2,
        true_rank_vector=true_rank_vector,
        allow_infeasible_diagnostic=allow_infeasible_diagnostic,
    )
