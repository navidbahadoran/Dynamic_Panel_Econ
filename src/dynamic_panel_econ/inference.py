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
from .lowrank import numerical_rank
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
    failure_code: str | None = None


@dataclass(slots=True)
class TangentBlock:
    u1: np.ndarray
    u0: np.ndarray
    v1: np.ndarray
    v0: np.ndarray
    rank: int
    n: int
    t: int


@dataclass(slots=True)
class RieszSystem:
    """Target-independent tangent coordinates and empirical Gram operator."""

    fitted: Coefficients
    design: Design
    ranks: tuple[int, ...]
    blocks: list[TangentBlock]
    offsets: np.ndarray
    operator: LinearOperator
    spectrum: dict[str, Any] = field(default_factory=dict)

    @property
    def dimension(self) -> int:
        return int(self.offsets[-1])

    def from_coordinates(self, vector: np.ndarray) -> Coefficients:
        matrices = []
        for index, block in enumerate(self.blocks):
            values = vector[self.offsets[index] : self.offsets[index + 1]]
            rank = block.rank
            k_size = rank * rank
            left_size = (block.n - rank) * rank
            core = values[:k_size].reshape(rank, rank)
            left = values[k_size : k_size + left_size].reshape(block.n - rank, rank)
            right = values[k_size + left_size :].reshape(rank, block.t - rank)
            matrices.append(
                block.u1 @ core @ block.v1.T
                + block.u0 @ left @ block.v1.T
                + block.u1 @ right @ block.v0.T
            )
        return from_matrices(matrices, len(self.fitted.A), len(self.fitted.B))

    def to_coordinates(self, theta: Coefficients) -> np.ndarray:
        values = []
        for matrix, block in zip(theta.matrices(), self.blocks, strict=True):
            values.extend(
                [
                    (block.u1.T @ matrix @ block.v1).ravel(),
                    (block.u0.T @ matrix @ block.v1).ravel(),
                    (block.u1.T @ matrix @ block.v0).ravel(),
                ]
            )
        return np.concatenate(values) if values else np.empty(0)

    def project(self, theta: Coefficients) -> Coefficients:
        return self.from_coordinates(self.to_coordinates(theta))


@dataclass(slots=True)
class SplitFitRecord:
    kind: str
    part: int
    rows: np.ndarray
    cols: np.ndarray
    y: np.ndarray
    design: Design
    fit: FitResult
    residuals: np.ndarray
    riesz_system: RieszSystem
    diagnostics: dict[str, Any]


@dataclass(slots=True)
class SplitFitBundle:
    time_parts: tuple[np.ndarray, np.ndarray]
    unit_parts: tuple[np.ndarray, np.ndarray]
    records: list[SplitFitRecord]
    coefficient_fit_count: int = 4


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


def prepare_riesz_system(
    fitted: Coefficients,
    design: Design,
    ranks: tuple[int, ...],
    *,
    compute_tangent_gram: bool = False,
    tangent_gram_tolerance: float = 1e-5,
    tangent_gram_max_iter: int = 500,
) -> RieszSystem:
    """Build one reusable nonredundant tangent-coordinate Gram system."""

    blocks: list[TangentBlock] = []
    sizes: list[int] = []
    for matrix, rank in zip(fitted.matrices(), ranks, strict=True):
        n, t = matrix.shape
        u, _, vt = np.linalg.svd(matrix, full_matrices=True)
        v = vt.T
        blocks.append(
            TangentBlock(u[:, :rank], u[:, rank:], v[:, :rank], v[:, rank:], rank, n, t)
        )
        sizes.append(rank * (n + t - rank))
    offsets = np.cumsum([0, *sizes])
    placeholder = LinearOperator(
        (int(offsets[-1]), int(offsets[-1])),
        matvec=lambda vector: vector,
        dtype=np.float64,
    )
    system = RieszSystem(fitted, design, ranks, blocks, offsets, placeholder)

    def apply(vector: np.ndarray) -> np.ndarray:
        tangent = system.from_coordinates(vector)
        return system.to_coordinates(adjoint(fitted_values(tangent, design), design))

    system.operator = LinearOperator(
        (system.dimension, system.dimension),
        matvec=apply,
        rmatvec=apply,
        dtype=np.float64,
    )
    if compute_tangent_gram:
        system.spectrum = _compute_tangent_gram_spectrum(
            system,
            tolerance=tangent_gram_tolerance,
            max_iter=tangent_gram_max_iter,
        )
    return system


def solve_riesz(
    direction: Coefficients,
    fitted: Coefficients,
    design: Design,
    ranks: tuple[int, ...],
    *,
    residuals: np.ndarray | None = None,
    tolerance: float = 1e-8,
    max_iter: int = 1000,
    system: RieszSystem | None = None,
) -> RieszResult:
    """Solve one target RHS on a reusable empirical tangent-Gram system."""

    system = system or prepare_riesz_system(fitted, design, ranks)
    rhs = system.to_coordinates(direction)
    operator = system.operator
    iteration_count = 0

    def callback(_: np.ndarray) -> None:
        nonlocal iteration_count
        iteration_count += 1

    solution, info = cg(operator, rhs, rtol=tolerance, atol=0.0, maxiter=max_iter, callback=callback)
    solver = "cg"

    def inaccurate(vector: np.ndarray, code: int) -> bool:
        if code != 0 or not np.all(np.isfinite(vector)):
            return True
        relative = np.linalg.norm(operator @ vector - rhs) / max(np.linalg.norm(rhs), 1e-15)
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
    q = system.from_coordinates(solution)
    weights = fitted_values(q, design)
    equation = operator @ solution - rhs
    equation_residual = float(np.linalg.norm(equation) / max(np.linalg.norm(rhs), 1e-15))
    aq = operator @ solution
    denominator = float(np.vdot(solution, solution))
    rayleigh = (
        float(np.vdot(solution, aq) / denominator)
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


def _compute_tangent_gram_spectrum(
    system: RieszSystem,
    *,
    tolerance: float = 1e-5,
    max_iter: int = 500,
) -> dict[str, Any]:
    """Estimate extreme eigenvalues of one prepared tangent Gram operator."""

    dimension = system.dimension
    if dimension == 0:
        return {
            "tangent_gram_coordinate_dimension": 0,
            "tangent_gram_smallest_eigenvalue": float("nan"),
            "tangent_gram_largest_eigenvalue": float("nan"),
            "tangent_gram_condition_number": float("nan"),
            "tangent_gram_eigensolver_converged": False,
            "tangent_gram_eigensolver_status": "zero_dimensional_tangent_space",
        }

    operator = system.operator
    if dimension == 1:
        smallest = largest = float((operator @ np.ones(1))[0])
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


def tangent_gram_spectrum(
    fitted: Coefficients,
    design: Design,
    ranks: tuple[int, ...],
    *,
    tolerance: float = 1e-5,
    max_iter: int = 500,
) -> dict[str, Any]:
    """Compatibility wrapper for one nonredundant tangent-Gram calculation."""

    system = prepare_riesz_system(fitted, design, ranks)
    return _compute_tangent_gram_spectrum(system, tolerance=tolerance, max_iter=max_iter)


def _empty_riesz(template: Coefficients, design: Design, target_norm: float) -> RieszResult:
    q = from_matrices(
        [np.zeros_like(matrix) for matrix in template.matrices()],
        len(template.A),
        len(template.B),
    )
    return RieszResult(
        q,
        np.zeros(template.shape),
        False,
        0,
        float("nan"),
        float("nan"),
        None,
        "not_called",
        target_norm,
    )


def _gram_failure_code(system: RieszSystem, minimum_eigenvalue_floor: float) -> str | None:
    if not system.spectrum:
        return None
    if not bool(system.spectrum["tangent_gram_eigensolver_converged"]):
        return "tangent_gram_eigensolver_failure"
    minimum = float(system.spectrum["tangent_gram_smallest_eigenvalue"])
    if not np.isfinite(minimum) or minimum < minimum_eigenvalue_floor:
        return "tangent_gram_nearly_singular"
    return None


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
    tangent_gram_min_eigenvalue_floor: float = 1e-10,
    target_support_tolerance: float = 1e-12,
    riesz_system: RieszSystem | None = None,
) -> InferenceResult:
    system = riesz_system or prepare_riesz_system(
        fit.theta,
        design,
        fit.ranks,
        compute_tangent_gram=compute_tangent_gram,
        tangent_gram_tolerance=tangent_gram_tolerance,
        tangent_gram_max_iter=tangent_gram_max_iter,
    )
    residuals = y - fitted_values(fit.theta, design)
    target_norm = float(np.linalg.norm(system.to_coordinates(direction)))
    cutoff = spatial_cutoff(*y.shape, c_sp) if spatial else 0
    diagnostics: dict[str, Any] = {
        "spatial_cutoff": cutoff,
        "target_tangent_norm": target_norm,
        "target_supported": bool(
            np.isfinite(target_norm) and target_norm > target_support_tolerance
        ),
        **system.spectrum,
    }
    if not diagnostics["target_supported"]:
        return InferenceResult(
            estimate=target_value(direction, fit.theta),
            standard_error=float("nan"),
            variance=float("nan"),
            riesz=_empty_riesz(fit.theta, design, target_norm),
            diagnostics=diagnostics,
            failure_code="target_unsupported_selected_rank",
        )
    gram_failure = _gram_failure_code(system, tangent_gram_min_eigenvalue_floor)
    if gram_failure is not None:
        return InferenceResult(
            estimate=target_value(direction, fit.theta),
            standard_error=float("nan"),
            variance=float("nan"),
            riesz=_empty_riesz(fit.theta, design, target_norm),
            diagnostics=diagnostics,
            failure_code=gram_failure,
        )
    riesz = solve_riesz(
        direction,
        fit.theta,
        design,
        fit.ranks,
        residuals=residuals,
        tolerance=riesz_tolerance,
        max_iter=riesz_max_iter,
        system=system,
    )
    scores = riesz.weights * residuals
    variance = bartlett_quadratic(scores, cutoff) if spatial else float(np.sum(scores * scores))
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


def _split_rank_diagnostics(
    fit: FitResult,
    relative_rank_floor: float,
) -> tuple[tuple[int, ...], list[dict[str, Any]], bool]:
    numerical_ranks = tuple(numerical_rank(matrix) for matrix in fit.theta.matrices())
    diagnostics = []
    rank_supported = True
    for label, matrix, supplied_rank in zip(
        ("A", "B", "H"), fit.theta.matrices(), fit.ranks, strict=True
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
        supported = supplied_rank == 0 or relative >= relative_rank_floor
        rank_supported = rank_supported and supported
        diagnostics.append(
            {
                "block": label,
                "supplied_rank": supplied_rank,
                "sigma_1": sigma_1,
                "sigma_r": sigma_r,
                "sigma_r_over_sigma_1": relative,
                "relative_rank_floor": relative_rank_floor,
                "computational_rank_supported": supported,
            }
        )
    return numerical_ranks, diagnostics, rank_supported


def prepare_split_fits(
    full_fit: FitResult,
    y: np.ndarray,
    design: Design,
    groups: np.ndarray,
    *,
    time_seed: int | np.random.SeedSequence,
    unit_seed: int | np.random.SeedSequence,
    split_relative_rank_floor: float = 1e-10,
    compute_tangent_gram: bool = False,
    tangent_gram_tolerance: float = 1e-5,
    tangent_gram_max_iter: int = 500,
    fit_options: dict[str, Any] | None = None,
) -> SplitFitBundle:
    """Compute exactly four target-independent split coefficient fits."""

    n, t = y.shape
    all_rows, all_cols = np.arange(n), np.arange(t)
    time_parts, unit_parts = balanced_partitions(
        n, t, groups, np.random.default_rng(time_seed), np.random.default_rng(unit_seed)
    )
    fit_options = fit_options or {}
    records: list[SplitFitRecord] = []

    def fit_one(rows: np.ndarray, cols: np.ndarray, kind: str, part: int) -> None:
        sub_design = subset_design(design, rows, cols)
        sub_y = y[np.ix_(rows, cols)]
        sub_initial = subset_coefficients(full_fit.theta, rows, cols)
        sub_fit = fit_fixed_rank(
            sub_y,
            sub_design,
            full_fit.ranks,
            initial=sub_initial,
            seed=part + (0 if kind == "time" else 2),
            **fit_options,
        )
        residual = sub_y - fitted_values(sub_fit.theta, sub_design)
        numerical_ranks, rank_diagnostics, rank_supported = _split_rank_diagnostics(
            sub_fit, split_relative_rank_floor
        )
        system = prepare_riesz_system(
            sub_fit.theta,
            sub_design,
            sub_fit.ranks,
            compute_tangent_gram=compute_tangent_gram,
            tangent_gram_tolerance=tangent_gram_tolerance,
            tangent_gram_max_iter=tangent_gram_max_iter,
        )
        diagnostics = {
                "kind": kind,
                "part": part,
                "n": len(rows),
                "t": len(cols),
                "converged": sub_fit.converged,
                "stationarity_residual": sub_fit.stationarity_residual,
                "max_envelope_ratio": sub_fit.max_envelope_ratio,
                "numerical_rank_vector": numerical_ranks,
                "rank_supported": rank_supported,
                "split_rank_singular_values": rank_diagnostics,
                **system.spectrum,
            }
        records.append(
            SplitFitRecord(
                kind,
                part,
                rows,
                cols,
                sub_y,
                sub_design,
                sub_fit,
                residual,
                system,
                diagnostics,
            )
        )

    for part, cols in enumerate(time_parts):
        fit_one(all_rows, cols, "time", part)
    for part, rows in enumerate(unit_parts):
        fit_one(rows, all_cols, "unit", part)
    return SplitFitBundle(time_parts, unit_parts, records, len(records))


def infer_corrected_target(
    direction: Coefficients,
    full_fit: FitResult,
    full_system: RieszSystem,
    y: np.ndarray,
    design: Design,
    split_bundle: SplitFitBundle,
    *,
    spatial: bool,
    c_sp: float = 1.0,
    riesz_tolerance: float = 1e-8,
    riesz_max_iter: int = 1000,
    target_rayleigh_floor: float = 1e-12,
    target_support_tolerance: float = 1e-12,
    tangent_gram_min_eigenvalue_floor: float = 1e-10,
) -> InferenceResult:
    """Infer one broad target using four already-computed split fits."""

    full = infer_target(
        direction,
        full_fit,
        y,
        design,
        spatial=spatial,
        c_sp=c_sp,
        riesz_tolerance=riesz_tolerance,
        riesz_max_iter=riesz_max_iter,
        target_support_tolerance=target_support_tolerance,
        tangent_gram_min_eigenvalue_floor=tangent_gram_min_eigenvalue_floor,
        riesz_system=full_system,
    )
    base_split_diagnostics = [dict(record.diagnostics) for record in split_bundle.records]
    common = {
        **full.diagnostics,
        "phi_full": full.estimate,
        "plugin_estimate": full.estimate,
        "plugin_standard_error": full.standard_error,
        "plugin_variance": full.variance,
        "split_coefficient_fit_count": split_bundle.coefficient_fit_count,
        "time_parts": [part.tolist() for part in split_bundle.time_parts],
        "unit_parts": [part.tolist() for part in split_bundle.unit_parts],
    }
    if full.failure_code is not None:
        full.corrected = True
        full.diagnostics.update({**common, "split_fits": base_split_diagnostics})
        return full

    weights_time = np.zeros_like(y)
    residuals_time = np.zeros_like(y)
    weights_unit = np.zeros_like(y)
    residuals_unit = np.zeros_like(y)
    phi_time = 0.0
    phi_unit = 0.0
    split_diagnostics: list[dict[str, Any]] = []
    for record in split_bundle.records:
        sub_direction = restrict_direction(direction, record.rows, record.cols)
        support_norm = float(np.linalg.norm(record.riesz_system.to_coordinates(sub_direction)))
        target_supported = bool(
            np.isfinite(support_norm) and support_norm > target_support_tolerance
        )
        diagnostics = {
            **record.diagnostics,
            "target_support_norm": support_norm,
            "target_supported": target_supported,
        }
        split_diagnostics.append(diagnostics)
        if not target_supported:
            return InferenceResult(
                estimate=float("nan"),
                standard_error=float("nan"),
                variance=float("nan"),
                riesz=full.riesz,
                corrected=True,
                diagnostics={**common, "split_fits": split_diagnostics},
                failure_code="split_target_unsupported_selected_rank",
            )
        gram_failure = _gram_failure_code(
            record.riesz_system, tangent_gram_min_eigenvalue_floor
        )
        if gram_failure is not None:
            return InferenceResult(
                estimate=float("nan"),
                standard_error=float("nan"),
                variance=float("nan"),
                riesz=full.riesz,
                corrected=True,
                diagnostics={**common, "split_fits": split_diagnostics},
                failure_code=f"split_{gram_failure}",
            )
        riesz = solve_riesz(
            sub_direction,
            record.fit.theta,
            record.design,
            record.fit.ranks,
            residuals=record.residuals,
            tolerance=riesz_tolerance,
            max_iter=riesz_max_iter,
            system=record.riesz_system,
        )
        diagnostics.update(
            {
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
        value = target_value(sub_direction, record.fit.theta)
        if record.kind == "time":
            phi_time += value
            weights_time[np.ix_(record.rows, record.cols)] = riesz.weights
            residuals_time[np.ix_(record.rows, record.cols)] = record.residuals
        else:
            phi_unit += value
            weights_unit[np.ix_(record.rows, record.cols)] = riesz.weights
            residuals_unit[np.ix_(record.rows, record.cols)] = record.residuals

    scores = corrected_scores(
        full.riesz.weights,
        y - fitted_values(full_fit.theta, design),
        weights_time,
        residuals_time,
        weights_unit,
        residuals_unit,
    )
    cutoff = spatial_cutoff(*y.shape, c_sp) if spatial else 0
    variance = bartlett_quadratic(scores, cutoff) if spatial else float(np.sum(scores**2))
    estimate = 3.0 * full.estimate - phi_time - phi_unit
    return InferenceResult(
        estimate=estimate,
        standard_error=float(np.sqrt(max(variance, 0.0))),
        variance=variance,
        riesz=full.riesz,
        corrected=True,
        diagnostics={
            **common,
            "spatial_cutoff": cutoff,
            "phi_time_sum": phi_time,
            "phi_unit_sum": phi_unit,
            "split_fits": split_diagnostics,
            "distinct_split_residual_scores": True,
        },
    )


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
    tangent_gram_min_eigenvalue_floor: float = 1e-10,
    fit_options: dict[str, Any] | None = None,
) -> InferenceResult:
    """Compatibility wrapper; new Monte Carlo code prepares the bundle once."""

    full_system = prepare_riesz_system(
        full_fit.theta,
        design,
        full_fit.ranks,
        compute_tangent_gram=compute_tangent_gram,
        tangent_gram_tolerance=tangent_gram_tolerance,
        tangent_gram_max_iter=tangent_gram_max_iter,
    )
    bundle = prepare_split_fits(
        full_fit,
        y,
        design,
        groups,
        time_seed=time_seed,
        unit_seed=unit_seed,
        split_relative_rank_floor=split_relative_rank_floor,
        compute_tangent_gram=compute_tangent_gram,
        tangent_gram_tolerance=tangent_gram_tolerance,
        tangent_gram_max_iter=tangent_gram_max_iter,
        fit_options=fit_options,
    )
    return infer_corrected_target(
        direction,
        full_fit,
        full_system,
        y,
        design,
        bundle,
        spatial=spatial,
        c_sp=c_sp,
        riesz_tolerance=riesz_tolerance,
        riesz_max_iter=riesz_max_iter,
        target_rayleigh_floor=target_rayleigh_floor,
        target_support_tolerance=target_support_tolerance,
        tangent_gram_min_eigenvalue_floor=tangent_gram_min_eigenvalue_floor,
    )
