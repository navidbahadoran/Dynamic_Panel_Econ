"""Joint supplied-rank least-squares and convex nuclear screening."""

from __future__ import annotations

import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import LinearConstraint, minimize, nnls

from .core import (
    Coefficients,
    Design,
    adjoint,
    fitted_values,
    frobenius_norm,
    from_matrices,
    max_abs,
    scale,
    zeros_for_design,
)
from .lowrank import singular_values, tangent_project, truncated_matrix

Array = NDArray[np.float64]
_FIXED_FIT_OBSERVER: ContextVar[Callable[[FitResult], None] | None] = ContextVar(
    "fixed_fit_observer", default=None
)
_NUCLEAR_FIT_OBSERVER: ContextVar[Callable[[NuclearFit], None] | None] = ContextVar(
    "nuclear_fit_observer", default=None
)


@dataclass(slots=True)
class FactorBlock:
    loading: Array
    factor: Array

    def matrix(self) -> Array:
        return self.loading @ self.factor.T


@dataclass(slots=True)
class FitResult:
    theta: Coefficients
    ranks: tuple[int, ...]
    objective: float
    converged: bool
    iterations: int
    objective_history: list[float]
    stationarity_residual: float
    max_envelope_ratio: float
    factors: list[FactorBlock] = field(repr=False)
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class NuclearFit:
    theta: Coefficients
    penalty: float
    objective: float
    converged: bool
    iterations: int
    objective_history: list[float]
    singular_values: list[list[float]]
    runtime_seconds: float | None = None


@contextmanager
def observe_fixed_rank_fits(
    observer: Callable[[FitResult], None],
) -> Iterator[None]:
    """Collect every actual supplied-rank numerical fit in the current worker."""

    token = _FIXED_FIT_OBSERVER.set(observer)
    try:
        yield
    finally:
        _FIXED_FIT_OBSERVER.reset(token)


@contextmanager
def observe_nuclear_fits(
    observer: Callable[[NuclearFit], None],
) -> Iterator[None]:
    token = _NUCLEAR_FIT_OBSERVER.set(observer)
    try:
        yield
    finally:
        _NUCLEAR_FIT_OBSERVER.reset(token)


def _to_theta(blocks: list[FactorBlock], p: int, k: int) -> Coefficients:
    return from_matrices([block.matrix() for block in blocks], p, k)


def _initial_blocks(
    shape: tuple[int, int],
    ranks: tuple[int, ...],
    rng: np.random.Generator,
    initial: Coefficients | None,
) -> list[FactorBlock]:
    n, t = shape
    matrices = initial.matrices() if initial is not None else [None] * len(ranks)
    blocks: list[FactorBlock] = []
    for rank, matrix in zip(ranks, matrices, strict=True):
        if rank == 0:
            blocks.append(FactorBlock(np.zeros((n, 0)), np.zeros((t, 0))))
            continue
        if matrix is None:
            loading = rng.normal(scale=0.1, size=(n, rank))
            factor = rng.normal(scale=0.1, size=(t, rank))
        else:
            u, singular, vt = np.linalg.svd(matrix, full_matrices=False)
            available = min(rank, singular.size)
            roots = np.sqrt(np.maximum(singular[:available], 0.0))
            loading = np.zeros((n, rank), dtype=np.float64)
            factor = np.zeros((t, rank), dtype=np.float64)
            loading[:, :available] = u[:, :available] * roots
            factor[:, :available] = vt[:available].T * roots
            if available < rank or np.linalg.matrix_rank(matrix) < rank:
                loading += rng.normal(scale=1e-4, size=loading.shape)
                factor += rng.normal(scale=1e-4, size=factor.shape)
        blocks.append(FactorBlock(loading.astype(np.float64), factor.astype(np.float64)))
    return blocks


def _renormalize(block: FactorBlock) -> None:
    if block.loading.shape[1] == 0:
        return
    q, r = np.linalg.qr(block.loading, mode="reduced")
    block.loading = q
    block.factor = block.factor @ r.T


def _stationarity(y: Array, design: Design, theta: Coefficients, ranks: tuple[int, ...]) -> float:
    n, t = y.shape
    residual = y - fitted_values(theta, design)
    gradient = scale(adjoint(-residual, design), 1.0 / (n * t))
    projected = tangent_project(gradient, theta, ranks)
    denominator = max(float(np.linalg.norm(y)) / np.sqrt(n * t), 1.0)
    return frobenius_norm(projected) / denominator


def _fit_fixed_rank_unconstrained(
    y: Array,
    design: Design,
    ranks: tuple[int, ...],
    *,
    initial: Coefficients | None = None,
    seed: int | np.random.SeedSequence = 0,
    max_sweeps: int = 200,
    objective_rtol: float = 1e-8,
    stationarity_tol: float = 1e-6,
    coefficient_bound: float = 10.0,
    lstsq_rcond: float = 1e-10,
    diagnostic_context: str | None = None,
) -> FitResult:
    """Joint ALS for all supplied-rank coefficient matrices."""

    started = time.perf_counter()
    n, t = y.shape
    if design.shape != y.shape or len(ranks) != len(design.y_lags) + len(design.x) + 1:
        raise ValueError("incompatible data, design, and rank vector")
    if any(rank < 0 or rank > min(n, t) for rank in ranks):
        raise ValueError("invalid supplied rank")
    rng = np.random.default_rng(seed)
    blocks = _initial_blocks((n, t), ranks, rng, initial)
    regressors = design.regressors()
    initial_envelope = max_abs(_to_theta(blocks, len(design.y_lags), len(design.x)))
    envelope_history = [initial_envelope]

    def objective() -> float:
        residual = y - fitted_values(_to_theta(blocks, len(design.y_lags), len(design.x)), design)
        return float(np.vdot(residual, residual) / (2.0 * n * t))

    history = [objective()]
    converged = False
    for _sweep in range(1, max_sweeps + 1):
        for i in range(n):
            pieces = [
                regressors[j][i, :, None] * block.factor
                for j, block in enumerate(blocks)
                if block.loading.shape[1]
            ]
            if pieces:
                solution, *_ = np.linalg.lstsq(np.concatenate(pieces, axis=1), y[i], rcond=lstsq_rcond)
                offset = 0
                for block in blocks:
                    rank = block.loading.shape[1]
                    if rank:
                        block.loading[i] = solution[offset : offset + rank]
                        offset += rank
        for block in blocks:
            _renormalize(block)
        for column in range(t):
            pieces = [
                regressors[j][:, column, None] * block.loading
                for j, block in enumerate(blocks)
                if block.factor.shape[1]
            ]
            if pieces:
                solution, *_ = np.linalg.lstsq(
                    np.concatenate(pieces, axis=1), y[:, column], rcond=lstsq_rcond
                )
                offset = 0
                for block in blocks:
                    rank = block.factor.shape[1]
                    if rank:
                        block.factor[column] = solution[offset : offset + rank]
                        offset += rank
        value = objective()
        envelope_history.append(
            max_abs(_to_theta(blocks, len(design.y_lags), len(design.x)))
        )
        if value > history[-1] + 1e-11 * max(1.0, history[-1]):
            break
        history.append(value)
        relative = abs(history[-2] - value) / max(abs(history[-2]), 1e-14)
        if relative <= objective_rtol:
            converged = True
            break
    theta = _to_theta(blocks, len(design.y_lags), len(design.x))
    stationarity = _stationarity(y, design, theta, ranks)
    converged = bool(converged or stationarity <= stationarity_tol)
    envelope_ratio = max_abs(theta) / coefficient_bound
    # Objective convergence and the stationarity diagnostic are stored separately.
    result = FitResult(
        theta=theta,
        ranks=ranks,
        objective=history[-1],
        converged=converged,
        iterations=len(history) - 1,
        objective_history=history,
        stationarity_residual=stationarity,
        max_envelope_ratio=envelope_ratio,
        factors=blocks,
        diagnostics={
            "stationarity_pass": stationarity <= stationarity_tol,
            "bound_active": envelope_ratio >= 1.0,
            "singular_values": singular_values(theta),
            "runtime_seconds": time.perf_counter() - started,
            "diagnostic_context": diagnostic_context,
            "initial_coefficient_envelope": initial_envelope,
            "final_coefficient_envelope": max_abs(theta),
            "coefficient_envelope_history": envelope_history,
        },
    )
    return result


def _box_subproblem_matrices(
    blocks: list[FactorBlock],
    regressors: tuple[Array, ...],
    index: int,
    *,
    update_loading: bool,
) -> tuple[Array, Array, Array]:
    ranks = [block.loading.shape[1] for block in blocks]
    total_rank = sum(ranks)
    if update_loading:
        design = np.concatenate(
            [
                regressors[j][index, :, None] * block.factor
                for j, block in enumerate(blocks)
                if ranks[j]
            ],
            axis=1,
        )
        current = np.concatenate(
            [block.loading[index] for block in blocks if block.loading.shape[1]]
        )
        constraint_rows = sum(block.factor.shape[0] for block in blocks if block.factor.shape[1])
    else:
        design = np.concatenate(
            [
                regressors[j][:, index, None] * block.loading
                for j, block in enumerate(blocks)
                if ranks[j]
            ],
            axis=1,
        )
        current = np.concatenate(
            [block.factor[index] for block in blocks if block.factor.shape[1]]
        )
        constraint_rows = sum(
            block.loading.shape[0] for block in blocks if block.loading.shape[1]
        )
    constraints = np.zeros((constraint_rows, total_rank), dtype=np.float64)
    row_offset = 0
    rank_offset = 0
    for block, rank in zip(blocks, ranks, strict=True):
        if not rank:
            continue
        basis = block.factor if update_loading else block.loading
        width = basis.shape[0]
        constraints[row_offset : row_offset + width, rank_offset : rank_offset + rank] = basis
        row_offset += width
        rank_offset += rank
    return design, current, constraints


def _linear_box_kkt_residual(
    design: Array,
    outcome: Array,
    constraints: Array,
    solution: Array,
    bound: float,
    active_tolerance: float,
) -> float:
    scale_factor = max(1, outcome.size)
    gradient = design.T @ (design @ solution - outcome) / scale_factor
    values = constraints @ solution
    upper = values >= bound - active_tolerance
    lower = values <= -bound + active_tolerance
    normals = np.concatenate((constraints[upper], -constraints[lower]), axis=0)
    if normals.size:
        multipliers, _ = nnls(normals.T, -gradient)
        lagrangian_gradient = gradient + normals.T @ multipliers
    else:
        lagrangian_gradient = gradient
    denominator = max(
        1.0,
        float(np.linalg.norm(design.T @ outcome / scale_factor)),
    )
    return float(np.linalg.norm(lagrangian_gradient) / denominator)


def _solve_linear_box_subproblem(
    design: Array,
    outcome: Array,
    constraints: Array,
    initial: Array,
    bound: float,
    *,
    max_iterations: int,
    tolerance: float,
    constraint_tolerance: float,
) -> tuple[Array, bool, str, int]:
    scale_factor = max(1, outcome.size)

    def objective(value: Array) -> float:
        residual = design @ value - outcome
        return float(np.vdot(residual, residual) / (2.0 * scale_factor))

    def gradient(value: Array) -> Array:
        return design.T @ (design @ value - outcome) / scale_factor

    initial_objective = objective(initial)
    result = minimize(
        objective,
        initial,
        jac=gradient,
        constraints=LinearConstraint(constraints, -bound, bound),
        method="SLSQP",
        options={"maxiter": max_iterations, "ftol": tolerance, "disp": False},
    )
    candidate = np.asarray(result.x, dtype=np.float64)
    violation = max(0.0, float(np.max(np.abs(constraints @ candidate))) - bound)
    kkt = _linear_box_kkt_residual(
        design,
        outcome,
        constraints,
        candidate,
        bound,
        max(constraint_tolerance * 10.0, 1e-8),
    )
    nonincreasing = objective(candidate) <= initial_objective + tolerance * max(
        1.0, abs(initial_objective)
    )
    accepted = bool(
        np.all(np.isfinite(candidate))
        and violation <= constraint_tolerance
        and nonincreasing
        and (result.success or kkt <= max(1e-6, np.sqrt(tolerance)))
    )
    return candidate, accepted, str(result.message), int(result.nit)


def _factor_kkt_residual(
    blocks: list[FactorBlock],
    regressors: tuple[Array, ...],
    y: Array,
    bound: float,
    active_tolerance: float,
) -> float:
    residuals = []
    for i in range(y.shape[0]):
        design, current, constraints = _box_subproblem_matrices(
            blocks, regressors, i, update_loading=True
        )
        residuals.append(
            _linear_box_kkt_residual(
                design, y[i], constraints, current, bound, active_tolerance
            )
        )
    for column in range(y.shape[1]):
        design, current, constraints = _box_subproblem_matrices(
            blocks, regressors, column, update_loading=False
        )
        residuals.append(
            _linear_box_kkt_residual(
                design, y[:, column], constraints, current, bound, active_tolerance
            )
        )
    return float(max(residuals, default=0.0))


def _fit_fixed_rank_constrained(
    y: Array,
    design: Design,
    ranks: tuple[int, ...],
    *,
    initial: Coefficients,
    seed: int | np.random.SeedSequence,
    max_sweeps: int,
    objective_rtol: float,
    coefficient_bound: float,
    constraint_tolerance: float,
    constrained_kkt_tolerance: float,
    constrained_subproblem_tolerance: float,
    constrained_subproblem_max_iterations: int,
    diagnostic_context: str | None,
) -> FitResult:
    started = time.perf_counter()
    n, t = y.shape
    rng = np.random.default_rng(seed)
    initial_max = max_abs(initial)
    if initial_max >= coefficient_bound:
        start = scale(
            initial,
            coefficient_bound * (1.0 - 10.0 * constraint_tolerance) / initial_max,
        )
    else:
        start = initial
    blocks = _initial_blocks((n, t), ranks, rng, start)
    regressors = design.regressors()

    def objective() -> float:
        theta = _to_theta(blocks, len(design.y_lags), len(design.x))
        residual = y - fitted_values(theta, design)
        return float(np.vdot(residual, residual) / (2.0 * n * t))

    history = [objective()]
    messages: list[str] = []
    subproblem_iterations = 0
    solver_failed = False
    kkt = float("inf")
    for _sweep in range(1, max_sweeps + 1):
        for i in range(n):
            sub_design, current, constraints = _box_subproblem_matrices(
                blocks, regressors, i, update_loading=True
            )
            x, accepted, message, iterations = _solve_linear_box_subproblem(
                sub_design,
                y[i],
                constraints,
                current,
                coefficient_bound,
                max_iterations=constrained_subproblem_max_iterations,
                tolerance=constrained_subproblem_tolerance,
                constraint_tolerance=constraint_tolerance,
            )
            if not accepted:
                messages.append(f"loading[{i}]: {message}")
                solver_failed = True
                break
            offset = 0
            for block in blocks:
                rank = block.loading.shape[1]
                if rank:
                    block.loading[i] = x[offset : offset + rank]
                    offset += rank
            subproblem_iterations += iterations
        if solver_failed:
            break
        for block in blocks:
            _renormalize(block)
        for column in range(t):
            sub_design, current, constraints = _box_subproblem_matrices(
                blocks, regressors, column, update_loading=False
            )
            x, accepted, message, iterations = _solve_linear_box_subproblem(
                sub_design,
                y[:, column],
                constraints,
                current,
                coefficient_bound,
                max_iterations=constrained_subproblem_max_iterations,
                tolerance=constrained_subproblem_tolerance,
                constraint_tolerance=constraint_tolerance,
            )
            if not accepted:
                messages.append(f"factor[{column}]: {message}")
                solver_failed = True
                break
            offset = 0
            for block in blocks:
                rank = block.factor.shape[1]
                if rank:
                    block.factor[column] = x[offset : offset + rank]
                    offset += rank
            subproblem_iterations += iterations
        value = objective()
        history.append(value)
        if solver_failed:
            break
        kkt = _factor_kkt_residual(
            blocks,
            regressors,
            y,
            coefficient_bound,
            max(10.0 * constraint_tolerance, 1e-8),
        )
        relative = abs(history[-2] - value) / max(1.0, abs(history[-2]))
        if relative <= objective_rtol and kkt <= constrained_kkt_tolerance:
            break
    theta = _to_theta(blocks, len(design.y_lags), len(design.x))
    maximum = max_abs(theta)
    violation = max(0.0, maximum - coefficient_bound)
    finite = bool(np.isfinite(history[-1]) and np.all([np.all(np.isfinite(m)) for m in theta.matrices()]))
    feasibility_pass = finite and violation <= constraint_tolerance
    kkt_pass = np.isfinite(kkt) and kkt <= constrained_kkt_tolerance
    converged = bool(not solver_failed and feasibility_pass and kkt_pass)
    if not finite:
        status = "nonfinite_constrained_solution"
    elif not feasibility_pass:
        status = "constrained_feasibility_failure"
    elif solver_failed:
        status = "constrained_solver_failure"
    elif not kkt_pass:
        status = "constrained_optimality_failure"
    else:
        status = "success"
    return FitResult(
        theta=theta,
        ranks=ranks,
        objective=history[-1],
        converged=converged,
        iterations=len(history) - 1,
        objective_history=history,
        stationarity_residual=kkt,
        max_envelope_ratio=maximum / coefficient_bound,
        factors=blocks,
        diagnostics={
            "stationarity_pass": kkt_pass,
            "stationarity_type": "factor_space_box_KKT",
            "bound_active": maximum >= coefficient_bound - 10.0 * constraint_tolerance,
            "boundary_active": maximum >= coefficient_bound - 10.0 * constraint_tolerance,
            "singular_values": singular_values(theta),
            "runtime_seconds": time.perf_counter() - started,
            "constrained_runtime": time.perf_counter() - started,
            "diagnostic_context": diagnostic_context,
            "initial_coefficient_envelope": initial_max,
            "final_coefficient_envelope": maximum,
            "coefficient_envelope_history": [
                max_abs(start),
                maximum,
            ],
            "max_constraint_violation": violation,
            "constraint_tolerance": constraint_tolerance,
            "constrained_KKT_residual": kkt,
            "constrained_iterations": subproblem_iterations,
            "constrained_solver_status": status,
            "constrained_objective": history[-1],
            "constrained_algorithm": "alternating_exact_linear_box_QP",
            "constrained_messages": messages,
        },
    )


def fit_fixed_rank(
    y: Array,
    design: Design,
    ranks: tuple[int, ...],
    *,
    initial: Coefficients | None = None,
    seed: int | np.random.SeedSequence = 0,
    max_sweeps: int = 200,
    objective_rtol: float = 1e-8,
    stationarity_tol: float = 1e-6,
    coefficient_bound: float = 10.0,
    lstsq_rcond: float = 1e-10,
    interior_numerical_tolerance: float = 1e-8,
    constraint_tolerance: float = 1e-8,
    constrained_kkt_tolerance: float = 1e-5,
    constrained_subproblem_tolerance: float = 1e-10,
    constrained_subproblem_max_iterations: int = 200,
    diagnostic_context: str | None = None,
) -> FitResult:
    """Literal fixed-rank box-constrained LS with an interior ALS fast path."""

    unconstrained = _fit_fixed_rank_unconstrained(
        y,
        design,
        ranks,
        initial=initial,
        seed=seed,
        max_sweeps=max_sweeps,
        objective_rtol=objective_rtol,
        stationarity_tol=stationarity_tol,
        coefficient_bound=coefficient_bound,
        lstsq_rcond=lstsq_rcond,
        diagnostic_context=diagnostic_context,
    )
    unconstrained_max = max_abs(unconstrained.theta)
    inside = unconstrained_max < coefficient_bound - interior_numerical_tolerance
    if inside:
        unconstrained.diagnostics.update(
            {
                "unconstrained_max_abs": unconstrained_max,
                "unconstrained_inside_box": True,
                "unconstrained_outside_box": False,
                "constrained_fallback_used": False,
                "boundary_active": False,
                "max_constraint_violation": 0.0,
                "constrained_KKT_residual": None,
                "constrained_iterations": 0,
                "constrained_runtime": 0.0,
                "constrained_solver_status": "not_needed_interior_fast_path",
                "constrained_objective": None,
            }
        )
        result = unconstrained
    else:
        result = _fit_fixed_rank_constrained(
            y,
            design,
            ranks,
            initial=unconstrained.theta,
            seed=seed,
            max_sweeps=max_sweeps,
            objective_rtol=objective_rtol,
            coefficient_bound=coefficient_bound,
            constraint_tolerance=constraint_tolerance,
            constrained_kkt_tolerance=constrained_kkt_tolerance,
            constrained_subproblem_tolerance=constrained_subproblem_tolerance,
            constrained_subproblem_max_iterations=constrained_subproblem_max_iterations,
            diagnostic_context=diagnostic_context,
        )
        result.diagnostics.update(
            {
                "unconstrained_objective": unconstrained.objective,
                "unconstrained_max_abs": unconstrained_max,
                "unconstrained_inside_box": False,
                "unconstrained_outside_box": True,
                "constrained_fallback_used": True,
            }
        )
    observer = _FIXED_FIT_OBSERVER.get()
    if observer is not None:
        observer(result)
    return result


def _svt(matrix: Array, threshold: float) -> Array:
    u, singular, vt = np.linalg.svd(matrix, full_matrices=False)
    keep = np.maximum(singular - threshold, 0.0)
    return (u * keep) @ vt


def _nuclear_box_prox(
    matrix: Array,
    threshold: float,
    bound: float,
    max_iter: int,
    tolerance: float,
) -> Array:
    """Dykstra prox for nuclear norm plus the entrywise box indicator."""

    current = matrix.copy()
    p = np.zeros_like(matrix)
    q = np.zeros_like(matrix)
    for _ in range(max_iter):
        previous = current
        nuclear = _svt(current + p, threshold)
        p = current + p - nuclear
        current = np.clip(nuclear + q, -bound, bound)
        q = nuclear + q - current
        if np.linalg.norm(current - previous) <= tolerance * max(1.0, np.linalg.norm(previous)):
            break
    return current


def penalty_weights(design: Design) -> tuple[float, ...]:
    weights = [float(np.sqrt(np.mean(z * z))) for z in [*design.y_lags, *design.x]]
    return (*weights, 1.0)


def lambda_maximum(y: Array, design: Design) -> float:
    n, t = y.shape
    gradient = scale(adjoint(-y, design), 1.0 / (n * t))
    weights = penalty_weights(design)
    values = [np.linalg.svd(g, compute_uv=False)[0] / w for g, w in zip(gradient.matrices(), weights, strict=True)]
    return float(max(values))


def _nuclear_objective(
    y: Array, design: Design, theta: Coefficients, penalty: float, weights: tuple[float, ...]
) -> float:
    n, t = y.shape
    residual = y - fitted_values(theta, design)
    loss = float(np.vdot(residual, residual) / (2.0 * n * t))
    nuclear = sum(w * np.linalg.svd(m, compute_uv=False).sum() for m, w in zip(theta.matrices(), weights, strict=True))
    return loss + penalty * float(nuclear)


def fit_nuclear(
    y: Array,
    design: Design,
    penalty: float,
    *,
    initial: Coefficients | None = None,
    coefficient_bound: float = 10.0,
    max_iter: int = 500,
    tolerance: float = 1e-7,
    dykstra_max_iter: int = 100,
    dykstra_tolerance: float = 1e-9,
) -> NuclearFit:
    """Monotone proximal-gradient solution of the convex screening problem."""

    started = time.perf_counter()
    n, t = y.shape
    theta = initial.copy() if initial is not None else zeros_for_design(design)
    weights = penalty_weights(design)
    lipschitz = max(1.0, sum(float(np.max(z * z)) for z in design.regressors()) / (n * t))
    history = [_nuclear_objective(y, design, theta, penalty, weights)]
    converged = False
    for _iteration in range(1, max_iter + 1):
        residual = fitted_values(theta, design) - y
        gradient = scale(adjoint(residual, design), 1.0 / (n * t))
        while True:
            proposed = []
            for matrix, grad, weight in zip(theta.matrices(), gradient.matrices(), weights, strict=True):
                proposed.append(
                    _nuclear_box_prox(
                        matrix - grad / lipschitz,
                        penalty * weight / lipschitz,
                        coefficient_bound,
                        dykstra_max_iter,
                        dykstra_tolerance,
                    )
                )
            candidate = from_matrices(proposed, len(theta.A), len(theta.B))
            value = _nuclear_objective(y, design, candidate, penalty, weights)
            if value <= history[-1] + 1e-12:
                break
            lipschitz *= 2.0
        history.append(value)
        difference = np.sqrt(sum(np.vdot(a - b, a - b) for a, b in zip(candidate.matrices(), theta.matrices(), strict=True)))
        theta = candidate
        if difference <= tolerance * max(1.0, frobenius_norm(theta)):
            converged = True
            break
        lipschitz = max(lipschitz * 0.9, 1e-12)
    result = NuclearFit(
        theta,
        penalty,
        history[-1],
        converged,
        _iteration,
        history,
        singular_values(theta),
        time.perf_counter() - started,
    )
    observer = _NUCLEAR_FIT_OBSERVER.get()
    if observer is not None:
        observer(result)
    return result


def nuclear_path(
    y: Array,
    design: Design,
    *,
    gamma: float = 0.8,
    epsilon: float = 0.01,
    **kwargs: Any,
) -> list[NuclearFit]:
    maximum = lambda_maximum(y, design)
    length = 1 + int(np.ceil(np.log(epsilon) / np.log(gamma)))
    fits: list[NuclearFit] = []
    initial = None
    for index in range(length):
        fit = fit_nuclear(y, design, maximum * gamma**index, initial=initial, **kwargs)
        fits.append(fit)
        initial = fit.theta
    return fits


def adapt_initial(theta: Coefficients, ranks: tuple[int, ...]) -> Coefficients:
    matrices = [truncated_matrix(matrix, rank) for matrix, rank in zip(theta.matrices(), ranks, strict=True)]
    return from_matrices(matrices, len(theta.A), len(theta.B))
