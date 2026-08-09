"""Four cumulative Monte Carlo data-generating processes."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .core import Coefficients, Design
from .seeds import rng_for

Array = NDArray[np.float64]
SQRT3 = float(np.sqrt(3.0))
INITIAL_CONDITIONS = {
    "y_i_minus_50": 0.0,
    "x_i_minus_50": 0.0,
    "g_a_minus_50": 0.0,
    "g_b_minus_50": 0.0,
    "g_h_minus_50": 0.0,
    "f_x_minus_50": 0.0,
}
DGP4_TRUTH_NAMES = (
    "A_G1_fixed_time_true",
    "A_G2_fixed_time_true",
    "A_G2_minus_G1_fixed_time_true",
    "B_G1_fixed_time_true",
    "B_G2_fixed_time_true",
    "B_G2_minus_G1_fixed_time_true",
    "A_G1_time_average_true",
    "A_G2_time_average_true",
    "A_G2_minus_G1_time_average_true",
    "B_G1_time_average_true",
    "B_G2_time_average_true",
    "B_G2_minus_G1_time_average_true",
    "A_G2_minus_G1_fixed_time_raw_true",
    "A_G2_minus_G1_time_average_raw_true",
)


@dataclass(frozen=True, slots=True)
class DGPParameters:
    burn_in: int = 50
    rho_g: float = 0.5
    rho_s: float = 0.5
    rho_x: float = 0.5
    rho_fx: float = 0.5
    delta_x: float = 0.5
    eta_x: float = 0.3
    mu_f_a: float = 0.5
    kappa_f_a: float = 0.1
    mu_f_b: float = 0.6
    kappa_f_b: float = 0.15
    mu_lambda_a_1: float = 0.9
    mu_lambda_a_2: float = 1.1
    sigma_lambda_a: float = 0.08
    mu_lambda_b_1: float = 0.8
    mu_lambda_b_2: float = 1.2
    sigma_lambda_b: float = 0.25
    stability_bound: float = 0.85

    @classmethod
    def from_mapping(cls, values: dict[str, Any]) -> DGPParameters:
        names = cls.__dataclass_fields__.keys()
        return cls(**{key: values[key] for key in names if key in values})


@dataclass(slots=True)
class PanelData:
    y: Array
    design: Design
    theta0: Coefficients
    u: Array
    u_tilde: Array
    u_tilde_lag: Array
    groups: NDArray[np.int_]
    diagnostics: dict[str, Any] = field(default_factory=dict)
    truths: dict[str, float] = field(default_factory=dict)


@dataclass(slots=True)
class _RawDraw:
    a: Array
    a_raw: Array
    beta: Array
    h_raw: Array
    x: Array
    u_tilde: Array
    groups: NDArray[np.int_]
    c_a: float
    pre_max_a: float
    loading_a: Array
    loading_b: Array
    theoretical_raw_envelopes: dict[str, float] = field(default_factory=dict)
    stress_rescale_factors: dict[str, float] = field(default_factory=dict)


def _bounded(rng: np.random.Generator, shape: int | tuple[int, ...]) -> Array:
    return rng.uniform(-SQRT3, SQRT3, size=shape).astype(np.float64)


def _ar_uniform(rng: np.random.Generator, length: int, rho: float) -> Array:
    values = np.empty(length, dtype=np.float64)
    previous = INITIAL_CONDITIONS["g_a_minus_50"]
    scale = np.sqrt(1.0 - rho * rho)
    for j in range(length):
        previous = rho * previous + scale * float(_bounded(rng, 1)[0])
        values[j] = previous
    return values


def bounded_ar_envelope(rho: float) -> float:
    """Infinite-horizon deterministic envelope for the bounded AR recursion."""

    if abs(rho) >= 1.0:
        raise ValueError("bounded AR envelope requires abs(rho) < 1")
    return float(SQRT3 * np.sqrt(1.0 - rho * rho) / (1.0 - abs(rho)))


def rank_one_raw_envelopes(dgp: int, params: DGPParameters) -> dict[str, float]:
    """Deterministic raw coefficient envelopes implied by bounded primitives."""

    g_bound = bounded_ar_envelope(params.rho_g)
    if dgp == 4:
        loading_a = max(abs(params.mu_lambda_a_1), abs(params.mu_lambda_a_2)) + (
            params.sigma_lambda_a * SQRT3
        )
        loading_b = max(abs(params.mu_lambda_b_1), abs(params.mu_lambda_b_2)) + (
            params.sigma_lambda_b * SQRT3
        )
    else:
        loading_a = 1.0 + 0.1 * SQRT3
        loading_b = 1.0 + 0.4 * SQRT3
    factor_a = abs(params.mu_f_a) + abs(params.kappa_f_a) * g_bound
    factor_b = abs(params.mu_f_b) + abs(params.kappa_f_b) * g_bound
    h_bound = SQRT3 * bounded_ar_envelope(params.rho_g)
    return {
        "A_raw": float(loading_a * factor_a),
        "A": float(min(params.stability_bound, loading_a * factor_a)),
        "B": float(loading_b * factor_b),
        "H_raw": float(h_bound),
    }


def stress_rescale_factor(
    base_envelope: float,
    rank: int,
    component_strengths: tuple[float, ...],
) -> float:
    """Common higher-rank rescaling that preserves the rank-one envelope."""

    if rank <= 1:
        return 1.0
    extra = sum(
        3.0 * abs(component_strengths[min(component - 1, len(component_strengths) - 1)])
        for component in range(1, rank)
    )
    return float(base_envelope / (base_envelope + extra))


def coefficient_envelopes(
    dgp: int,
    params: DGPParameters,
    *,
    c_h: float,
    c_xi: float,
    ranks: tuple[int, int, int] = (1, 1, 1),
) -> dict[str, float]:
    """Deterministic final coefficient envelopes at fixed calibration constants."""

    raw = rank_one_raw_envelopes(dgp, params)
    values = {
        "A": raw["A"] if ranks[0] else 0.0,
        "B": raw["B"] if ranks[1] else 0.0,
        "H": float(abs(c_h * c_xi) * raw["H_raw"]) if ranks[2] else 0.0,
    }
    values["all"] = max(values.values())
    return values


def _spatial_errors(
    rng: np.random.Generator, n: int, length: int, rho_s: float, spatial: bool
) -> tuple[Array, Array]:
    sigma2 = rng.uniform(0.5, 1.5, size=n)
    sigma = np.sqrt(sigma2)
    eps = rng.normal(size=(n, length)).astype(np.float64)
    if not spatial:
        return sigma[:, None] * eps, sigma2
    z = np.empty_like(eps)
    z[0] = eps[0]
    innovation_scale = np.sqrt(1.0 - rho_s * rho_s)
    for i in range(1, n):
        z[i] = rho_s * z[i - 1] + innovation_scale * eps[i]
    return sigma[:, None] * z, sigma2


def _draw_raw(
    dgp: int,
    n: int,
    t: int,
    rng: np.random.Generator,
    params: DGPParameters,
) -> _RawDraw:
    if dgp not in {1, 2, 3, 4}:
        raise ValueError("dgp must be 1, 2, 3, or 4")
    if dgp == 4 and n % 2:
        raise ValueError("DGP 4 requires even N for equal deterministic groups")
    length = t + params.burn_in
    groups = np.zeros(n, dtype=np.int64)
    groups[n // 2 :] = 1

    z_a, z_b = _bounded(rng, n), _bounded(rng, n)
    if dgp == 4:
        lambda_a = np.where(groups == 0, params.mu_lambda_a_1, params.mu_lambda_a_2)
        lambda_b = np.where(groups == 0, params.mu_lambda_b_1, params.mu_lambda_b_2)
        lambda_a = lambda_a + params.sigma_lambda_a * z_a
        lambda_b = lambda_b + params.sigma_lambda_b * z_b
    else:
        lambda_a = 1.0 + 0.1 * z_a
        lambda_b = 1.0 + 0.4 * z_b

    g_a = _ar_uniform(rng, length, params.rho_g)
    g_b = _ar_uniform(rng, length, params.rho_g)
    f_a = params.mu_f_a + params.kappa_f_a * g_a
    f_b = params.mu_f_b + params.kappa_f_b * g_b
    a_raw = lambda_a[:, None] * f_a[None, :]
    pre_max = float(np.max(np.abs(a_raw)))
    c_a = min(1.0, params.stability_bound / pre_max)
    a = c_a * a_raw
    beta = lambda_b[:, None] * f_b[None, :]

    lambda_h = _bounded(rng, n)
    g_h = _ar_uniform(rng, length, params.rho_g)
    h_raw = lambda_h[:, None] * g_h[None, :]

    # Include the primitive disturbance at t=-50 for predetermined x at t=-49.
    u_all, _ = _spatial_errors(rng, n, length + 1, params.rho_s, dgp >= 2)
    lambda_x = _bounded(rng, n)
    sigma_e2 = rng.uniform(0.5, 1.5, size=n)
    e = rng.normal(size=(n, length)) * np.sqrt(sigma_e2)[:, None]
    f_x = _ar_uniform(rng, length, params.rho_fx)
    x = np.empty((n, length), dtype=np.float64)
    previous_x = np.full(n, INITIAL_CONDITIONS["x_i_minus_50"], dtype=np.float64)
    x_innovation_scale = np.sqrt(1.0 - params.rho_x * params.rho_x)
    for j in range(length):
        predetermined = params.eta_x * u_all[:, j] if dgp >= 3 else 0.0
        current = (
            params.rho_x * previous_x
            + params.delta_x * lambda_x * f_x[j]
            + predetermined
            + x_innovation_scale * e[:, j]
        )
        x[:, j] = current
        previous_x = current
    return _RawDraw(
        a=a,
        a_raw=a_raw,
        beta=beta,
        h_raw=h_raw,
        x=x,
        u_tilde=u_all[:, 1:],
        groups=groups,
        c_a=c_a,
        pre_max_a=pre_max,
        loading_a=lambda_a,
        loading_b=lambda_b,
        theoretical_raw_envelopes=rank_one_raw_envelopes(dgp, params),
        stress_rescale_factors={"A": 1.0, "B": 1.0, "H": 1.0},
    )


def _draw_rank_stress_raw(
    dgp: int,
    n: int,
    t: int,
    ranks: tuple[int, int, int],
    rng: np.random.Generator,
    params: DGPParameters,
    component_strengths: tuple[float, ...],
) -> _RawDraw:
    """Draw exact-rank stress matrices rescaled to rank-one support envelopes."""

    raw = _draw_raw(dgp, n, t, rng, params)
    length = t + params.burn_in
    base_envelopes = rank_one_raw_envelopes(dgp, params)

    def expand(base: Array, rank: int, label: str) -> tuple[Array, float]:
        if rank == 0:
            return np.zeros_like(base), 1.0
        matrix = base.copy()
        for component in range(1, rank):
            strength = component_strengths[min(component - 1, len(component_strengths) - 1)]
            loading = _bounded(rng, n)
            factor = _bounded(rng, length)
            matrix += float(strength) * loading[:, None] * factor[None, :]
        factor = stress_rescale_factor(base_envelopes[label], rank, component_strengths)
        return factor * matrix, factor

    a_raw, scale_a = expand(raw.a_raw, ranks[0], "A_raw")
    beta, scale_b = expand(raw.beta, ranks[1], "B")
    h_raw, scale_h = expand(raw.h_raw, ranks[2], "H_raw")
    pre_max = float(np.max(np.abs(a_raw))) if ranks[0] else 0.0
    c_a = min(1.0, params.stability_bound / pre_max) if pre_max else 1.0
    return replace(
        raw,
        a=c_a * a_raw,
        a_raw=a_raw,
        beta=beta,
        h_raw=h_raw,
        c_a=c_a,
        pre_max_a=pre_max,
        theoretical_raw_envelopes=base_envelopes,
        stress_rescale_factors={"A": scale_a, "B": scale_b, "H": scale_h},
    )


def _simulate(raw: _RawDraw, c_h: float, c_xi: float) -> tuple[Array, Array, Array]:
    n, length = raw.x.shape
    h0 = c_xi * c_h * raw.h_raw
    u = c_xi * raw.u_tilde
    y = np.empty((n, length), dtype=np.float64)
    previous = np.full(n, INITIAL_CONDITIONS["y_i_minus_50"], dtype=np.float64)
    for j in range(length):
        current = raw.a[:, j] * previous + raw.beta[:, j] * raw.x[:, j] + h0[:, j] + u[:, j]
        y[:, j] = current
        previous = current
    return y, h0, u


def group_truths(theta: Coefficients, a_raw: Array, groups: NDArray[np.int_]) -> dict[str, float]:
    """Exact realized DGP-4 group truths after scaling and raw A contrasts."""

    n, t = theta.shape
    fixed_index = max(1, t // 2) - 1
    g1, g2 = groups == 0, groups == 1
    out: dict[str, float] = {}
    for label, matrix in (("A", theta.A[0]), ("B", theta.B[0])):
        fixed1 = float(matrix[g1, fixed_index].mean())
        fixed2 = float(matrix[g2, fixed_index].mean())
        avg1 = float(matrix[g1].mean())
        avg2 = float(matrix[g2].mean())
        out[f"{label}_G1_fixed_time_true"] = fixed1
        out[f"{label}_G2_fixed_time_true"] = fixed2
        out[f"{label}_G2_minus_G1_fixed_time_true"] = fixed2 - fixed1
        out[f"{label}_G1_time_average_true"] = avg1
        out[f"{label}_G2_time_average_true"] = avg2
        out[f"{label}_G2_minus_G1_time_average_true"] = avg2 - avg1
    raw1 = float(a_raw[g1, fixed_index].mean())
    raw2 = float(a_raw[g2, fixed_index].mean())
    out["A_G2_minus_G1_fixed_time_raw_true"] = raw2 - raw1
    out["A_G2_minus_G1_time_average_raw_true"] = float(a_raw[g2].mean() - a_raw[g1].mean())
    return out


def generate_panel(
    dgp: int,
    n: int,
    t: int,
    seed: int | np.random.SeedSequence,
    *,
    c_h: float = 1.0,
    c_xi: float = 1.0,
    params: DGPParameters | None = None,
    coefficient_bound: float | None = None,
    simulation_interior_margin: float = 1.0,
) -> PanelData:
    """Generate one observed panel after the prescribed 50-period burn-in."""

    params = params or DGPParameters()
    rng = np.random.default_rng(seed)
    raw = _draw_raw(dgp, n, t, rng, params)
    y_all, h_all, u_all = _simulate(raw, c_h, c_xi)
    start = params.burn_in
    observed = slice(start, start + t)
    lagged = slice(start - 1, start + t - 1)
    y = y_all[:, observed]
    design = Design([y_all[:, lagged]], [raw.x[:, observed]])
    theta0 = Coefficients([raw.a[:, observed]], [raw.beta[:, observed]], h_all[:, observed])
    envelopes = coefficient_envelopes(dgp, params, c_h=c_h, c_xi=c_xi)
    realized = [float(np.max(np.abs(matrix))) for matrix in theta0.matrices()]
    diagnostics: dict[str, Any] = {
        "dgp": dgp,
        "c_a": raw.c_a,
        "c_H": float(c_h),
        "c_xi": float(c_xi),
        "pre_scaling_max_abs_a": raw.pre_max_a,
        "final_max_abs_a": float(np.max(np.abs(theta0.A[0]))),
        "mean_a": float(theta0.A[0].mean()),
        "mean_b": float(theta0.B[0].mean()),
        "min_a": float(theta0.A[0].min()),
        "max_a": float(theta0.A[0].max()),
        "min_b": float(theta0.B[0].min()),
        "max_b": float(theta0.B[0].max()),
        "mu_lambda_a_1": params.mu_lambda_a_1,
        "mu_lambda_a_2": params.mu_lambda_a_2,
        "mu_lambda_b_1": params.mu_lambda_b_1,
        "mu_lambda_b_2": params.mu_lambda_b_2,
        "theoretical_max_abs_A": envelopes["A"],
        "theoretical_max_abs_B": envelopes["B"],
        "theoretical_max_abs_H": envelopes["H"],
        "theoretical_coefficient_envelope": envelopes["all"],
        "realized_max_abs_A": realized[0],
        "realized_max_abs_B": realized[1],
        "realized_max_abs_H": realized[2],
        "realized_coefficient_envelope": max(realized),
        "coefficient_bound_B": coefficient_bound,
        "required_simulation_interior_margin": simulation_interior_margin,
        "deterministic_interior_margin": (
            float(coefficient_bound - envelopes["all"])
            if coefficient_bound is not None
            else float("nan")
        ),
        "coefficient_envelope_condition_pass": (
            bool(envelopes["all"] <= coefficient_bound - simulation_interior_margin)
            if coefficient_bound is not None
            else None
        ),
    }
    truths = group_truths(theta0, raw.a_raw[:, observed], raw.groups) if dgp == 4 else {}
    return PanelData(
        y=y,
        design=design,
        theta0=theta0,
        u=u_all[:, observed],
        u_tilde=raw.u_tilde[:, observed],
        u_tilde_lag=raw.u_tilde[:, slice(start - 1, start + t - 1)],
        groups=raw.groups,
        diagnostics=diagnostics,
        truths=truths,
    )


def generate_rank_stress_panel(
    dgp: int,
    n: int,
    t: int,
    ranks: tuple[int, int, int],
    seed: int | np.random.SeedSequence,
    *,
    component_strengths: tuple[float, ...] = (1.0, 1.0),
    c_h: float = 1.0,
    c_xi: float = 1.0,
    params: DGPParameters | None = None,
    coefficient_bound: float | None = None,
    simulation_interior_margin: float = 1.0,
) -> PanelData:
    """Generate the separate rank-selector stress design with exact configured ranks."""

    params = params or DGPParameters()
    rng = np.random.default_rng(seed)
    stress_raw = _draw_rank_stress_raw(
        dgp,
        n,
        t,
        ranks,
        rng,
        params,
        component_strengths,
    )
    y_all, h_all, u_all = _simulate(stress_raw, c_h, c_xi)
    start = params.burn_in
    observed = slice(start, start + t)
    lagged = slice(start - 1, start + t - 1)
    theta0 = Coefficients(
        [stress_raw.a[:, observed]], [stress_raw.beta[:, observed]], h_all[:, observed]
    )
    envelopes = coefficient_envelopes(
        dgp, params, c_h=c_h, c_xi=c_xi, ranks=ranks
    )
    realized = [float(np.max(np.abs(matrix))) for matrix in theta0.matrices()]
    return PanelData(
        y=y_all[:, observed],
        design=Design([y_all[:, lagged]], [stress_raw.x[:, observed]]),
        theta0=theta0,
        u=u_all[:, observed],
        u_tilde=stress_raw.u_tilde[:, observed],
        u_tilde_lag=stress_raw.u_tilde[:, slice(start - 1, start + t - 1)],
        groups=stress_raw.groups,
        diagnostics={
            "dgp": dgp,
            "true_rank_vector": ranks,
            "component_strengths": component_strengths,
            "c_a": stress_raw.c_a,
            "pre_scaling_max_abs_a": stress_raw.pre_max_a,
            "final_max_abs_a": float(np.max(np.abs(theta0.A[0]))),
            "stress_rescale_factor_A": stress_raw.stress_rescale_factors["A"],
            "stress_rescale_factor_B": stress_raw.stress_rescale_factors["B"],
            "stress_rescale_factor_H": stress_raw.stress_rescale_factors["H"],
            "theoretical_max_abs_A": envelopes["A"],
            "theoretical_max_abs_B": envelopes["B"],
            "theoretical_max_abs_H": envelopes["H"],
            "theoretical_coefficient_envelope": envelopes["all"],
            "realized_max_abs_A": realized[0],
            "realized_max_abs_B": realized[1],
            "realized_max_abs_H": realized[2],
            "realized_coefficient_envelope": max(realized),
            "coefficient_bound_B": coefficient_bound,
            "required_simulation_interior_margin": simulation_interior_margin,
            "deterministic_interior_margin": (
                float(coefficient_bound - envelopes["all"])
                if coefficient_bound is not None
                else float("nan")
            ),
            "coefficient_envelope_condition_pass": (
                bool(envelopes["all"] <= coefficient_bound - simulation_interior_margin)
                if coefficient_bound is not None
                else None
            ),
        },
    )


def group_gap_pilot(
    n: int,
    t: int,
    master_seed: int,
    params: DGPParameters,
    gap_grid: list[float],
    floor: float,
    auto_adjust: bool,
    draws: int = 20,
) -> tuple[DGPParameters, dict[str, Any]]:
    """Prespecification-only DGP-4 pilot; never conditions on production results."""

    center = 0.5 * (params.mu_lambda_a_1 + params.mu_lambda_a_2)
    original_gap = params.mu_lambda_a_2 - params.mu_lambda_a_1
    candidates = sorted(set([original_gap, *gap_grid])) if auto_adjust else [original_gap]
    reports: list[dict[str, float]] = []
    chosen = candidates[-1]
    for gap in candidates:
        trial = replace(
            params,
            mu_lambda_a_1=center - gap / 2.0,
            mu_lambda_a_2=center + gap / 2.0,
        )
        fixed1, fixed2, fixed, average1, average2, averaged, scaling = [], [], [], [], [], [], []
        for rep in range(draws):
            panel = generate_panel(
                4,
                n,
                t,
                rng_for(master_seed, "dgp4_group_gap_pilot", n, t, rep),
                params=trial,
            )
            fixed1.append(panel.truths["A_G1_fixed_time_true"])
            fixed2.append(panel.truths["A_G2_fixed_time_true"])
            fixed.append(panel.truths["A_G2_minus_G1_fixed_time_true"])
            average1.append(panel.truths["A_G1_time_average_true"])
            average2.append(panel.truths["A_G2_time_average_true"])
            averaged.append(panel.truths["A_G2_minus_G1_time_average_true"])
            scaling.append(panel.diagnostics["c_a"])
        report = {
            "gap": float(gap),
            "mean_G1_fixed_time_true": float(np.mean(fixed1)),
            "mean_G2_fixed_time_true": float(np.mean(fixed2)),
            "mean_fixed_time_difference": float(np.mean(fixed)),
            "mean_G1_time_average_true": float(np.mean(average1)),
            "mean_G2_time_average_true": float(np.mean(average2)),
            "mean_time_average_difference": float(np.mean(averaged)),
            "min_abs_mean_difference": float(min(abs(np.mean(fixed)), abs(np.mean(averaged)))),
            "mean_c_a": float(np.mean(scaling)),
        }
        reports.append(report)
        if report["min_abs_mean_difference"] >= floor:
            chosen = gap
            break
    resolved = replace(
        params,
        mu_lambda_a_1=center - chosen / 2.0,
        mu_lambda_a_2=center + chosen / 2.0,
    )
    return resolved, {
        "auto_adjust_enabled": bool(auto_adjust),
        "configured_floor": float(floor),
        "original_means": [params.mu_lambda_a_1, params.mu_lambda_a_2],
        "chosen_gap": float(chosen),
        "frozen_means": [resolved.mu_lambda_a_1, resolved.mu_lambda_a_2],
        "draws": int(draws),
        "candidate_reports": reports,
        "parameter_snapshot": asdict(resolved),
    }
