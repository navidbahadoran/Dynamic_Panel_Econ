"""Matrix-free empirical Riesz inference and exact split-panel correction."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, cg, eigsh, lsmr, minres

from .core import (
    Coefficients,
    Design,
    adjoint,
    fitted_values,
    from_matrices,
    subset_coefficients,
    subset_design,
)
from .estimation import FitResult, fit_fixed_rank
from .lowrank import numerical_rank, tangent_project
from .spatial import bartlett_quadratic, spatial_cutoff
from .targets import restrict_direction, target_value


@dataclass(slots=True)
class RieszResult:
    q: Coefficients
    weights: np.ndarray
    converged: bool
    iterations: int
    equation_residual: float
    target_rayleigh_quotient: float
    weighted_residual_identity: float | None = None
    solver: str = "cg"
    target_tangent_norm: float = float("nan")


@dataclass(slots=True)
class InferenceResult:
    estimate: float
    standard_error: float
    variance: float
    riesz: RieszResult
    corrected: bool = False
    diagnostics: dict[str, Any] = field(default_factory=dict)


def corrected_scores(
    full_weights: np.ndarray,
    full_residuals: np.ndarray,
    time_weights: np.ndarray,
    time_residuals: np.ndarray,
    unit_weights: np.ndarray,
    unit_residuals: np.ndarray,
) -> np.ndarray:
    """Observation score using three genuinely distinct fitted residual arrays."""

    return (
        3.0 * full_weights * full_residuals
        - time_weights * time_residuals
        - unit_weights * unit_residuals
    )


def _layout(template: Coefficients) -> tuple[list[tuple[int, int]], list[int]]:
    shapes = [matrix.shape for matrix in template.matrices()]
    sizes = [int(np.prod(shape)) for shape in shapes]
    return shapes, sizes


def _flatten(theta: Coefficients) -> np.ndarray:
    return np.concatenate([matrix.ravel() for matrix in theta.matrices()])


def _unflatten(vector: np.ndarray, template: Coefficients) -> Coefficients:
    shapes, sizes = _layout(template)
    splits = np.cumsum(sizes)[:-1]
    matrices = [part.reshape(shape) for part, shape in zip(np.split(vector, splits), shapes, strict=True)]
    return from_matrices(matrices, len(template.A), len(template.B))


def solve_riesz(
    direction: Coefficients,
    fitted: Coefficients,
    design: Design,
    ranks: tuple[int, ...],
    *,
    residuals: np.ndarray | None = None,
    tolerance: float = 1e-8,
    max_iter: int = 1000,
) -> RieszResult:
    """Solve the empirical tangent-space Riesz equation without a dense basis."""

    rhs_theta = tangent_project(direction, fitted, ranks)
    rhs = _flatten(rhs_theta)
    size = rhs.size

    def apply(vector: np.ndarray) -> np.ndarray:
        candidate = tangent_project(_unflatten(vector, fitted), fitted, ranks)
        image = fitted_values(candidate, design)
        normal = tangent_project(adjoint(image, design), fitted, ranks)
        return _flatten(normal)

    operator = LinearOperator((size, size), matvec=apply, rmatvec=apply, dtype=np.float64)
    iteration_count = 0

    def callback(_: np.ndarray) -> None:
        nonlocal iteration_count
        iteration_count += 1

    solution, info = cg(operator, rhs, rtol=tolerance, atol=0.0, maxiter=max_iter, callback=callback)
    solver = "cg"

    def inaccurate(vector: np.ndarray, code: int) -> bool:
        if code != 0 or not np.all(np.isfinite(vector)):
            return True
        relative = np.linalg.norm(apply(vector) - rhs) / max(np.linalg.norm(rhs), 1e-15)
        return bool(relative > max(10 * tolerance, 1e-7))

    if inaccurate(solution, info):
        iteration_count = 0
        solution, info = minres(operator, rhs, rtol=tolerance, maxiter=max_iter, callback=callback)
        solver = "minres"
    if inaccurate(solution, info):
        least_squares = lsmr(operator, rhs, atol=tolerance, btol=tolerance, maxiter=max_iter)
        solution = least_squares[0]
        info = 0 if least_squares[1] in {1, 2} else least_squares[1]
        iteration_count = int(least_squares[2])
        solver = "lsmr"
    q = tangent_project(_unflatten(solution, fitted), fitted, ranks)
    weights = fitted_values(q, design)
    equation = apply(_flatten(q)) - rhs
    equation_residual = float(np.linalg.norm(equation) / max(np.linalg.norm(rhs), 1e-15))
    projected_solution = _flatten(q)
    aq = apply(projected_solution)
    denominator = float(np.vdot(projected_solution, projected_solution))
    rayleigh = (
        float(np.vdot(projected_solution, aq) / denominator)
        if denominator > 0.0
        else float("nan")
    )
    identity = None
    if residuals is not None:
        identity = float(
            abs(np.sum(weights * residuals))
            / max(np.linalg.norm(weights) * np.linalg.norm(residuals), 1e-15)
        )
    converged = bool(
        np.isfinite(equation_residual)
        and equation_residual <= max(10 * tolerance, 1e-7)
    )
    return RieszResult(
        q,
        weights,
        converged,
        iteration_count,
        equation_residual,
        rayleigh,
        identity,
        solver,
        float(np.linalg.norm(rhs)),
    )


def tangent_gram_spectrum(
    fitted: Coefficients,
    design: Design,
    ranks: tuple[int, ...],
    *,
    tolerance: float = 1e-5,
    max_iter: int = 500,
) -> dict[str, Any]:
    """Estimate extreme eigenvalues on nonredundant tangent coordinates.

    The coordinate basis contains exactly ``r(n+t-r)`` orthonormal directions
    per rank-r block, so ambient normal-space zero eigenvalues are excluded.
    """

    blocks: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int]] = []
    sizes: list[int] = []
    for matrix, rank in zip(fitted.matrices(), ranks, strict=True):
        n, t = matrix.shape
        u, _, vt = np.linalg.svd(matrix, full_matrices=True)
        v = vt.T
        u1, u0 = u[:, :rank], u[:, rank:]
        v1, v0 = v[:, :rank], v[:, rank:]
        blocks.append((u1, u0, v1, v0, n, t))
        sizes.append(rank * (n + t - rank))
    dimension = int(sum(sizes))
    if dimension == 0:
        return {
            "tangent_gram_coordinate_dimension": 0,
            "tangent_gram_smallest_eigenvalue": float("nan"),
            "tangent_gram_largest_eigenvalue": float("nan"),
            "tangent_gram_condition_number": float("nan"),
            "tangent_gram_eigensolver_converged": False,
            "tangent_gram_eigensolver_status": "zero_dimensional_tangent_space",
        }

    offsets = np.cumsum([0, *sizes])

    def from_coordinates(vector: np.ndarray) -> Coefficients:
        matrices = []
        for index, (u1, u0, v1, v0, n, t) in enumerate(blocks):
            rank = ranks[index]
            values = vector[offsets[index] : offsets[index + 1]]
            k_size = rank * rank
            l_size = (n - rank) * rank
            k = values[:k_size].reshape(rank, rank)
            left_normal = values[k_size : k_size + l_size].reshape(n - rank, rank)
            m = values[k_size + l_size :].reshape(rank, t - rank)
            matrices.append(
                u1 @ k @ v1.T + u0 @ left_normal @ v1.T + u1 @ m @ v0.T
            )
        return from_matrices(matrices, len(fitted.A), len(fitted.B))

    def to_coordinates(theta: Coefficients) -> np.ndarray:
        values = []
        for matrix, (u1, u0, v1, v0, _, _) in zip(
            theta.matrices(), blocks, strict=True
        ):
            values.extend(
                [
                    (u1.T @ matrix @ v1).ravel(),
                    (u0.T @ matrix @ v1).ravel(),
                    (u1.T @ matrix @ v0).ravel(),
                ]
            )
        return np.concatenate(values)

    def apply(vector: np.ndarray) -> np.ndarray:
        tangent = from_coordinates(vector)
        return to_coordinates(adjoint(fitted_values(tangent, design), design))

    operator = LinearOperator(
        (dimension, dimension), matvec=apply, rmatvec=apply, dtype=np.float64
    )
    if dimension == 1:
        smallest = largest = float(apply(np.ones(1))[0])
        converged, status = True, "exact_one_dimensional"
    else:
        initial = np.full(dimension, 1.0 / np.sqrt(dimension))
        try:
            smallest = float(
                eigsh(
                    operator,
                    k=1,
                    which="SA",
                    return_eigenvectors=False,
                    tol=tolerance,
                    maxiter=max_iter,
                    v0=initial,
                )[0]
            )
            largest = float(
                eigsh(
                    operator,
                    k=1,
                    which="LA",
                    return_eigenvectors=False,
                    tol=tolerance,
                    maxiter=max_iter,
                    v0=initial,
                )[0]
            )
            converged, status = True, "converged"
        except ArpackNoConvergence as exc:
            values = np.asarray(exc.eigenvalues)
            smallest = float(np.min(values)) if values.size else float("nan")
            largest = float(np.max(values)) if values.size else float("nan")
            converged, status = False, "arpack_no_convergence"
    condition = (
        float(largest / smallest)
        if np.isfinite(smallest) and np.isfinite(largest) and smallest > 0.0
        else float("inf")
    )
    return {
        "tangent_gram_coordinate_dimension": dimension,
        "tangent_gram_smallest_eigenvalue": smallest,
        "tangent_gram_largest_eigenvalue": largest,
        "tangent_gram_condition_number": condition,
        "tangent_gram_eigensolver_converged": converged,
        "tangent_gram_eigensolver_status": status,
    }
def infer_target(
    direction: Coefficients,
    fit: FitResult,
    y: np.ndarray,
    design: Design,
    *,
    spatial: bool,
    c_sp: float = 1.0,
    riesz_tolerance: float = 1e-8,
    riesz_max_iter: int = 1000,
    compute_tangent_gram: bool = False,
    tangent_gram_tolerance: float = 1e-5,
    tangent_gram_max_iter: int = 500,
) -> InferenceResult:
    residuals = y - fitted_values(fit.theta, design)
    riesz = solve_riesz(
        direction,
        fit.theta,
        design,
        fit.ranks,
        residuals=residuals,
        tolerance=riesz_tolerance,
        max_iter=riesz_max_iter,
    )
    scores = riesz.weights * residuals
    cutoff = spatial_cutoff(*y.shape, c_sp) if spatial else 0
    variance = bartlett_quadratic(scores, cutoff) if spatial else float(np.sum(scores * scores))
    diagnostics: dict[str, Any] = {"spatial_cutoff": cutoff}
    if compute_tangent_gram:
        diagnostics.update(
            tangent_gram_spectrum(
                fit.theta,
                design,
                fit.ranks,
                tolerance=tangent_gram_tolerance,
                max_iter=tangent_gram_max_iter,
            )
        )
    return InferenceResult(
        estimate=target_value(direction, fit.theta),
        standard_error=float(np.sqrt(max(variance, 0.0))),
        variance=variance,
        riesz=riesz,
        diagnostics=diagnostics,
    )


def balanced_partitions(
    n: int,
    t: int,
    groups: np.ndarray,
    rng_time: np.random.Generator,
    rng_unit: np.random.Generator,
) -> tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]:
    times = rng_time.permutation(t)
    time_parts = tuple(np.sort(part) for part in np.array_split(times, 2))
    units: list[list[int]] = [[], []]
    for group in np.unique(groups):
        shuffled = rng_unit.permutation(np.flatnonzero(groups == group))
        pieces = np.array_split(shuffled, 2)
        units[0].extend(int(value) for value in pieces[0])
        units[1].extend(int(value) for value in pieces[1])
    unit_parts = (np.array(sorted(units[0]), dtype=int), np.array(sorted(units[1]), dtype=int))
    return (time_parts[0], time_parts[1]), unit_parts


def split_correct_target(
    direction: Coefficients,
    full_fit: FitResult,
    y: np.ndarray,
    design: Design,
    groups: np.ndarray,
    *,
    time_seed: int | np.random.SeedSequence,
    unit_seed: int | np.random.SeedSequence,
    spatial: bool,
    c_sp: float = 1.0,
    riesz_tolerance: float = 1e-8,
    riesz_max_iter: int = 1000,
    target_rayleigh_floor: float = 1e-12,
    target_support_tolerance: float = 1e-12,
    split_relative_rank_floor: float = 1e-10,
    compute_tangent_gram: bool = False,
    tangent_gram_tolerance: float = 1e-5,
    tangent_gram_max_iter: int = 500,
    fit_options: dict[str, Any] | None = None,
) -> InferenceResult:
    """Apply the prescribed four-fit two-way correction for a broad target."""

    n, t = y.shape
    all_rows, all_cols = np.arange(n), np.arange(t)
    time_parts, unit_parts = balanced_partitions(
        n, t, groups, np.random.default_rng(time_seed), np.random.default_rng(unit_seed)
    )
    fit_options = fit_options or {}
    weights_time = np.zeros_like(y)
    residuals_time = np.zeros_like(y)
    weights_unit = np.zeros_like(y)
    residuals_unit = np.zeros_like(y)
    phi_time = 0.0
    phi_unit = 0.0
    split_diagnostics: list[dict[str, Any]] = []

    def process(
        rows: np.ndarray, cols: np.ndarray, kind: str, part: int
    ) -> tuple[float, np.ndarray, np.ndarray, FitResult, RieszResult]:
        sub_design = subset_design(design, rows, cols)
        sub_y = y[np.ix_(rows, cols)]
        sub_initial = subset_coefficients(full_fit.theta, rows, cols)
        sub_direction = restrict_direction(direction, rows, cols)
        sub_fit = fit_fixed_rank(
            sub_y,
            sub_design,
            full_fit.ranks,
            initial=sub_initial,
            seed=part,
            **fit_options,
        )
        residual = sub_y - fitted_values(sub_fit.theta, sub_design)
        numerical_ranks = tuple(numerical_rank(matrix) for matrix in sub_fit.theta.matrices())
        split_rank_diagnostics = []
        rank_supported = True
        for label, matrix, supplied_rank in zip(
            ("A", "B", "H"), sub_fit.theta.matrices(), full_fit.ranks, strict=True
        ):
            singular_values = np.linalg.svd(matrix, compute_uv=False)
            sigma_1 = float(singular_values[0]) if singular_values.size else 0.0
            sigma_r = (
                float(singular_values[supplied_rank - 1])
                if supplied_rank > 0 and singular_values.size >= supplied_rank
                else 0.0
            )
            relative = (
                sigma_r / max(sigma_1, np.finfo(float).tiny)
                if supplied_rank > 0
                else 1.0
            )
            supported = supplied_rank == 0 or relative >= split_relative_rank_floor
            rank_supported = rank_supported and supported
            split_rank_diagnostics.append(
                {
                    "block": label,
                    "supplied_rank": supplied_rank,
                    "sigma_1": sigma_1,
                    "sigma_r": sigma_r,
                    "sigma_r_over_sigma_1": relative,
                    "relative_rank_floor": split_relative_rank_floor,
                    "computational_rank_supported": supported,
                }
            )
        projected_direction = tangent_project(
            sub_direction, sub_fit.theta, sub_fit.ranks
        )
        support_norm = float(
            np.sqrt(sum(np.vdot(matrix, matrix) for matrix in projected_direction.matrices()))
        )
        target_supported = bool(np.isfinite(support_norm) and support_norm > target_support_tolerance)
        riesz = solve_riesz(
            sub_direction,
            sub_fit.theta,
            sub_design,
            sub_fit.ranks,
            residuals=residual,
            tolerance=riesz_tolerance,
            max_iter=riesz_max_iter,
        )
        split_diagnostics.append(
            {
                "kind": kind,
                "part": part,
                "n": len(rows),
                "t": len(cols),
                "converged": sub_fit.converged,
                "stationarity_residual": sub_fit.stationarity_residual,
                "max_envelope_ratio": sub_fit.max_envelope_ratio,
                "numerical_rank_vector": numerical_ranks,
                "rank_supported": rank_supported,
                "split_rank_singular_values": split_rank_diagnostics,
                "target_support_norm": support_norm,
                "target_supported": target_supported,
                "riesz_converged": riesz.converged,
                "riesz_equation_residual": riesz.equation_residual,
                "riesz_target_rayleigh_quotient": riesz.target_rayleigh_quotient,
                "riesz_target_stable": bool(
                    np.isfinite(riesz.target_rayleigh_quotient)
                    and riesz.target_rayleigh_quotient >= target_rayleigh_floor
                ),
                "riesz_solver": riesz.solver,
            }
        )
        return target_value(sub_direction, sub_fit.theta), riesz.weights, residual, sub_fit, riesz

    for part, cols in enumerate(time_parts):
        value, weights, residuals, _, _ = process(all_rows, cols, "time", part)
        phi_time += value
        weights_time[np.ix_(all_rows, cols)] = weights
        residuals_time[np.ix_(all_rows, cols)] = residuals
    for part, rows in enumerate(unit_parts):
        value, weights, residuals, _, _ = process(rows, all_cols, "unit", part)
        phi_unit += value
        weights_unit[np.ix_(rows, all_cols)] = weights
        residuals_unit[np.ix_(rows, all_cols)] = residuals

    full = infer_target(
        direction,
        full_fit,
        y,
        design,
        spatial=spatial,
        c_sp=c_sp,
        riesz_tolerance=riesz_tolerance,
        riesz_max_iter=riesz_max_iter,
        compute_tangent_gram=compute_tangent_gram,
        tangent_gram_tolerance=tangent_gram_tolerance,
        tangent_gram_max_iter=tangent_gram_max_iter,
    )
    scores = corrected_scores(
        full.riesz.weights,
        y - fitted_values(full_fit.theta, design),
        weights_time,
        residuals_time,
        weights_unit,
        residuals_unit,
    )
    cutoff = spatial_cutoff(n, t, c_sp) if spatial else 0
    variance = bartlett_quadratic(scores, cutoff) if spatial else float(np.sum(scores**2))
    estimate = 3.0 * full.estimate - phi_time - phi_unit
    return InferenceResult(
        estimate=estimate,
        standard_error=float(np.sqrt(max(variance, 0.0))),
        variance=variance,
        riesz=full.riesz,
        corrected=True,
        diagnostics={
            **full.diagnostics,
            "spatial_cutoff": cutoff,
            "phi_full": full.estimate,
            "phi_time_sum": phi_time,
            "phi_unit_sum": phi_unit,
            "plugin_estimate": full.estimate,
            "plugin_standard_error": full.standard_error,
            "plugin_variance": full.variance,
            "time_parts": [part.tolist() for part in time_parts],
            "unit_parts": [part.tolist() for part in unit_parts],
            "split_fits": split_diagnostics,
            "distinct_split_residual_scores": True,
        },
    )
