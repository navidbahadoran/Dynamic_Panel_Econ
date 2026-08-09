"""Deterministic, production-separated DGP calibration."""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq

from .dgp import (
    DGPParameters,
    _draw_rank_stress_raw,
    _draw_raw,
    _RawDraw,
    coefficient_envelopes,
    rank_one_raw_envelopes,
    stress_rescale_factor,
)
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


def population_u_tilde_variance(dgp: int, params: DGPParameters) -> float:
    """Marginal population variance of the primitive idiosyncratic disturbance."""

    if dgp not in {1, 2, 3, 4}:
        raise ValueError("dgp must be 1, 2, 3, or 4")
    if abs(params.rho_s) >= 1.0:
        raise ValueError("spatial AR requires abs(rho_s) < 1")
    # E[sigma_i^2] = 1 for sigma_i^2 ~ U[0.5, 1.5].  For DGPs 2--4,
    # z_0 ~ N(0,1) and z_i=rho_s*z_{i-1}+sqrt(1-rho_s^2)*eps_i,
    # so every z_i has marginal variance one.  DGP 1 uses eps_i directly.
    return 1.0


def population_h_raw_variance(
    t: int,
    params: DGPParameters,
    *,
    rank: int = 1,
    component_strengths: tuple[float, ...] = (1.0, 1.0),
) -> float:
    """Average observed-entry population variance of the raw H matrix."""

    if t <= 0:
        raise ValueError("t must be positive")
    if rank <= 0:
        return 0.0
    rho = params.rho_g
    columns = np.arange(params.burn_in + 1, params.burn_in + t + 1, dtype=np.float64)
    base_variance = float(np.mean(1.0 - rho ** (2.0 * columns)))
    added_variance = sum(
        float(component_strengths[min(j - 1, len(component_strengths) - 1)]) ** 2
        for j in range(1, rank)
    )
    raw_envelope = rank_one_raw_envelopes(1, params)["H_raw"]
    rescale = stress_rescale_factor(raw_envelope, rank, component_strengths)
    return float(rescale * rescale * (base_variance + added_variance))


def deterministic_c_h(
    dgp: int,
    t: int,
    params: DGPParameters,
    *,
    pi_h: float,
    rank: int = 1,
    component_strengths: tuple[float, ...] = (1.0, 1.0),
) -> float:
    """Ex-ante H scale computed only from analytical population moments."""

    if not 0.0 < pi_h < 1.0:
        raise ValueError("pi_h must lie strictly between zero and one")
    var_u = population_u_tilde_variance(dgp, params)
    var_h = population_h_raw_variance(
        t, params, rank=rank, component_strengths=component_strengths
    )
    if var_h <= 0.0:
        raise ValueError("positive H rank is required for pi_h calibration")
    return float(np.sqrt((pi_h / (1.0 - pi_h)) * var_u / var_h))


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
    component_strengths: tuple[float, ...] = (1.0, 1.0),
    allow_infeasible_diagnostic: bool = False,
    tolerance: float = 1e-10,
) -> CalibrationResult:
    """Calibrate a fixed collection of baseline or actual stress-design raw draws."""

    draws = len(raws)
    start = params.burn_in
    observed = slice(start, start + t)
    var_u = float(np.mean([np.var(raw.u_tilde[:, observed]) for raw in raws]))
    var_h = float(np.mean([np.var(raw.h_raw[:, observed]) for raw in raws]))
    if var_u <= 0.0 or var_h <= 0.0:
        raise RuntimeError("nonpositive calibration variance")
    population_var_u = population_u_tilde_variance(dgp, params)
    population_var_h = population_h_raw_variance(
        t,
        params,
        rank=true_rank_vector[2],
        component_strengths=component_strengths,
    )
    c_h = deterministic_c_h(
        dgp,
        t,
        params,
        pi_h=pi_h,
        rank=true_rank_vector[2],
        component_strengths=component_strengths,
    )

    # Conditional on each raw draw, y(c_xi)=y_base+c_xi*y_scale.  Compute both
    # recursions once so bracketing does not rerun the dynamic simulation.
    r2_quadratics: list[tuple[float, float, float, float]] = []
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
        observed_base = y_base[:, observed]
        observed_scale = y_scale[:, observed]
        centered_base = observed_base - observed_base.mean()
        centered_scale = observed_scale - observed_scale.mean()
        primitive_u = raw.u_tilde[:, observed]
        r2_quadratics.append(
            (
                float(np.vdot(centered_base, centered_base)),
                float(2.0 * np.vdot(centered_base, centered_scale)),
                float(np.vdot(centered_scale, centered_scale)),
                float(np.vdot(primitive_u, primitive_u)),
            )
        )

    def average_r2(c_xi: float) -> float:
        values = [
            1.0 - c_xi * c_xi * u_ss / (base_ss + c_xi * cross + c_xi * c_xi * scale_ss)
            for base_ss, cross, scale_ss, u_ss in r2_quadratics
        ]
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
            float(brentq(root, *bracket, xtol=tolerance, rtol=tolerance))
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
    realized_calibration_h_share = (c_h * c_h * var_h) / (
        c_h * c_h * var_h + var_u
    )
    achieved_h_share = (c_h * c_h * population_var_h) / (
        c_h * c_h * population_var_h + population_var_u
    )
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
        "population_var_u_tilde": population_var_u,
        "population_var_h_raw": population_var_h,
        "population_pi_h": achieved_h_share,
        "realized_calibration_h_share": float(realized_calibration_h_share),
        "c_h_source": "analytical_population_moments",
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
    tolerance: float = 1e-10,
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
        tolerance=tolerance,
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
    tolerance: float = 1e-10,
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
        component_strengths=component_strengths,
        allow_infeasible_diagnostic=allow_infeasible_diagnostic,
        tolerance=tolerance,
    )


def load_frozen_calibrations(
    path: str | Path,
) -> dict[tuple[int, int, int, tuple[int, ...] | None], CalibrationResult]:
    """Load and validate an ex-ante frozen DGP-calibration table."""

    source = Path(path)
    with source.open("rb") as handle:
        payload = tomllib.load(handle)
    if int(payload.get("schema_version", 0)) != 1:
        raise ValueError("unsupported frozen calibration schema")
    results: dict[tuple[int, int, int, tuple[int, ...] | None], CalibrationResult] = {}
    for row in payload.get("calibration", []):
        is_stress = bool(row["rank_stress"])
        ranks = tuple(int(value) for value in row["true_rank_vector"])
        key = (int(row["dgp"]), int(row["n"]), int(row["t"]), ranks if is_stress else None)
        if key in results:
            raise ValueError(f"duplicate frozen calibration cell: {key}")
        numeric = {
            name: float(row[name])
            for name in (
                "c_h",
                "c_xi",
                "pi_h",
                "achieved_r2",
                "C_A",
                "C_beta",
                "C_H",
                "C_Theta",
            )
        }
        if not all(np.isfinite(value) for value in numeric.values()):
            raise ValueError(f"nonfinite frozen calibration value: {key}")
        if ranks[1] == 0 and numeric["c_xi"] != 1.0:
            raise ValueError(f"zero-slope rank requires c_xi=1: {key}")
        if not np.isclose(
            numeric["C_Theta"],
            max(numeric["C_A"], numeric["C_beta"], numeric["C_H"]),
            rtol=0.0,
            atol=1e-12,
        ):
            raise ValueError(f"inconsistent coefficient envelope: {key}")
        diagnostics = dict(row)
        diagnostics.update(
            {
                "theoretical_max_abs_A": numeric["C_A"],
                "theoretical_max_abs_B": numeric["C_beta"],
                "theoretical_max_abs_H": numeric["C_H"],
                "theoretical_coefficient_envelope": numeric["C_Theta"],
                "calibration_source": "frozen_ex_ante_table",
                "r2_scale_identified": ranks[1] > 0,
            }
        )
        results[key] = CalibrationResult(
            dgp=key[0],
            n=key[1],
            t=key[2],
            c_h=numeric["c_h"],
            c_xi=numeric["c_xi"],
            target_h_share=numeric["pi_h"],
            achieved_h_share=numeric["pi_h"],
            target_r2=(float(row["intended_r2"]) if ranks[1] > 0 else None),
            achieved_r2=numeric["achieved_r2"],
            diagnostics=diagnostics,
        )
    if not results:
        raise ValueError("frozen calibration table is empty")
    return results
