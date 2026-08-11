"""Semantically exact numerical solver for the Revision-10 cap+1 pilot."""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.sparse.linalg import LinearOperator, lsmr

from .core import Coefficients, Design, fitted_values, max_abs, scale
from .estimation import (
    _FIXED_FIT_OBSERVER,
    FactorBlock,
    FitResult,
    _box_subproblem_matrices,
    _factor_kkt_residual,
    _initial_blocks,
    _linear_box_kkt_residual,
    _renormalize,
    _stationarity,
    _to_theta,
)
from .lowrank import singular_values

Array = NDArray[np.float64]


@dataclass(slots=True)
class ActiveSetResult:
    solution: Array
    accepted: bool
    status: str
    iterations: int
    active_set: tuple[int, ...]
    kkt_residual: float
    condition_number: float


def balanced_factor_block(matrix: Array, width: int) -> FactorBlock:
    """Return a deterministic width-``width`` factorization of the same product.

    Exact zero singular directions have a zero loading and an orthonormal time-side
    direction.  The zero product is therefore preserved while the at-most-rank width
    remains available to a subsequent loading update.
    """

    if width < 0 or width > min(matrix.shape):
        raise ValueError("invalid factor width")
    if width == 0:
        return FactorBlock(
            np.zeros((matrix.shape[0], 0), dtype=np.float64),
            np.zeros((matrix.shape[1], 0), dtype=np.float64),
        )
    u, values, vt = np.linalg.svd(matrix, full_matrices=False)
    u = u[:, :width]
    values = values[:width]
    vt = vt[:width]
    roots = np.sqrt(values)
    loading = u * roots
    factor = vt.T * roots
    zero = values == 0.0
    if np.any(zero):
        loading[:, zero] = 0.0
        factor[:, zero] = vt[zero].T
    return FactorBlock(loading.astype(np.float64), factor.astype(np.float64))


def rebalance_blocks(blocks: list[FactorBlock]) -> list[FactorBlock]:
    """Balance every block without changing any coefficient product."""

    return [balanced_factor_block(block.matrix(), block.loading.shape[1]) for block in blocks]


def _copy_blocks(blocks: list[FactorBlock]) -> list[FactorBlock]:
    return [FactorBlock(block.loading.copy(), block.factor.copy()) for block in blocks]


def _objective(y: Array, design: Design, blocks: list[FactorBlock]) -> float:
    theta = _to_theta(blocks, len(design.y_lags), len(design.x))
    residual = y - fitted_values(theta, design)
    return float(np.vdot(residual, residual) / (2.0 * y.size))


def _condition_from_singular(values: Array) -> float:
    if values.size == 0:
        return 1.0
    positive = values[values > 0.0]
    if positive.size != values.size:
        return float("inf")
    return float(positive[0] / positive[-1])


def _stable_lstsq(design: Array, outcome: Array, rcond: float) -> tuple[Array, float]:
    """Solve the same least-squares problem in reversibly equilibrated coordinates."""

    norms = np.linalg.norm(design, axis=0)
    scales = np.ones_like(norms)
    nonzero = norms > 0.0
    scales[nonzero] = 1.0 / norms[nonzero]
    equilibrated = design * scales
    solution, _, _, singular = np.linalg.lstsq(equilibrated, outcome, rcond=rcond)
    return scales * solution, _condition_from_singular(singular)


def _active_set_qp(
    design: Array,
    outcome: Array,
    constraints: Array,
    initial: Array,
    bound: float,
    *,
    max_iterations: int,
    tolerance: float,
    constraint_tolerance: float,
    warm_active_set: tuple[int, ...] = (),
) -> ActiveSetResult:
    """Solve the existing convex quadratic/linear-box subproblem by active set."""

    scale_factor = max(1, outcome.size)
    norms = np.sqrt(
        np.sum(design * design, axis=0) / scale_factor
        + np.sum(constraints * constraints, axis=0) / max(1, constraints.shape[0])
    )
    coordinate_scale = np.ones_like(norms)
    nonzero = norms > 0.0
    coordinate_scale[nonzero] = 1.0 / norms[nonzero]
    x_design = design * coordinate_scale
    x_constraints = constraints * coordinate_scale
    z = initial / coordinate_scale
    inequalities = np.concatenate((x_constraints, -x_constraints), axis=0)
    bounds = np.full(inequalities.shape[0], bound, dtype=np.float64)
    values = inequalities @ z
    if np.max(values - bounds, initial=0.0) > constraint_tolerance:
        return ActiveSetResult(
            initial,
            False,
            "infeasible_initial_point",
            0,
            (),
            float("inf"),
            float("inf"),
        )
    hessian = x_design.T @ x_design / scale_factor
    linear = -(x_design.T @ outcome) / scale_factor
    hessian_singular = np.linalg.svd(hessian, compute_uv=False)
    condition = _condition_from_singular(hessian_singular)
    active = sorted(
        set(int(index) for index in warm_active_set if 0 <= index < len(bounds))
        | set(np.flatnonzero(bounds - values <= max(10.0 * constraint_tolerance, 1e-10)))
    )
    step_tolerance = max(tolerance, np.finfo(float).eps)
    status = "maximum_iterations"
    iterations = 0
    for _iteration in range(1, max_iterations + 1):
        iterations = _iteration
        gradient = hessian @ z + linear
        active_matrix = inequalities[active] if active else np.empty((0, z.size))
        if active:
            kkt_matrix = np.block(
                [
                    [hessian, active_matrix.T],
                    [active_matrix, np.zeros((len(active), len(active)))],
                ]
            )
            rhs = np.concatenate((-gradient, np.zeros(len(active))))
            answer, *_ = np.linalg.lstsq(kkt_matrix, rhs, rcond=1e-12)
            direction = answer[: z.size]
            multipliers = answer[z.size :]
        else:
            direction, *_ = np.linalg.lstsq(hessian, -gradient, rcond=1e-12)
            multipliers = np.empty(0, dtype=np.float64)
        if np.linalg.norm(direction) <= step_tolerance * max(1.0, np.linalg.norm(z)):
            if not active or np.min(multipliers, initial=0.0) >= -np.sqrt(step_tolerance):
                status = "success"
                break
            remove_position = int(np.argmin(multipliers))
            del active[remove_position]
            continue
        products = inequalities @ direction
        alpha = 1.0
        blocking: int | None = None
        active_set = set(active)
        for index in range(len(bounds)):
            if index in active_set or products[index] <= 0.0:
                continue
            candidate = (bounds[index] - values[index]) / products[index]
            if candidate < alpha:
                alpha = max(0.0, float(candidate))
                blocking = index
        z = z + alpha * direction
        values = inequalities @ z
        if blocking is not None and alpha < 1.0 - step_tolerance:
            active.append(blocking)
            active.sort()
    solution = coordinate_scale * z
    violation = max(0.0, float(np.max(np.abs(constraints @ solution))) - bound)
    kkt = _linear_box_kkt_residual(
        design,
        outcome,
        constraints,
        solution,
        bound,
        max(10.0 * constraint_tolerance, 1e-8),
    )
    initial_residual = design @ initial - outcome
    final_residual = design @ solution - outcome
    initial_objective = float(np.vdot(initial_residual, initial_residual) / (2.0 * scale_factor))
    final_objective = float(np.vdot(final_residual, final_residual) / (2.0 * scale_factor))
    accepted = bool(
        np.all(np.isfinite(solution))
        and violation <= constraint_tolerance
        and final_objective <= initial_objective + tolerance * max(1.0, abs(initial_objective))
        and (status == "success" or kkt <= max(1e-6, np.sqrt(tolerance)))
    )
    final_inequality_values = np.concatenate(
        (constraints @ solution, -(constraints @ solution))
    )
    original_active = tuple(
        int(index)
        for index in np.flatnonzero(
            bound - final_inequality_values <= max(10.0 * constraint_tolerance, 1e-8)
        )
    )
    return ActiveSetResult(
        solution,
        accepted,
        status,
        iterations,
        original_active,
        kkt,
        condition,
    )


def _pack_directions(blocks: list[FactorBlock]) -> tuple[list[tuple[int, int, int, int]], int]:
    layout: list[tuple[int, int, int, int]] = []
    offset = 0
    for block in blocks:
        loading_size = block.loading.size
        factor_size = block.factor.size
        layout.append((offset, offset + loading_size, offset + loading_size + factor_size, block.loading.shape[1]))
        offset += loading_size + factor_size
    return layout, offset


def _gauss_newton_direction(
    y: Array,
    design: Design,
    blocks: list[FactorBlock],
) -> list[FactorBlock]:
    regressors = design.regressors()
    layout, size = _pack_directions(blocks)
    n, t = y.shape

    def matvec(vector: Array) -> Array:
        fitted = np.zeros((n, t), dtype=np.float64)
        for regressor, block, (start, middle, end, width) in zip(
            regressors, blocks, layout, strict=True
        ):
            du = vector[start:middle].reshape(n, width)
            dv = vector[middle:end].reshape(t, width)
            fitted += regressor * (du @ block.factor.T + block.loading @ dv.T)
        return fitted.ravel()

    def rmatvec(vector: Array) -> Array:
        matrix = vector.reshape(n, t)
        result = np.empty(size, dtype=np.float64)
        for regressor, block, (start, middle, end, _width) in zip(
            regressors, blocks, layout, strict=True
        ):
            weighted = regressor * matrix
            result[start:middle] = (weighted @ block.factor).ravel()
            result[middle:end] = (weighted.T @ block.loading).ravel()
        return result

    operator = LinearOperator((n * t, size), matvec=matvec, rmatvec=rmatvec, dtype=np.float64)
    residual = (y - fitted_values(_to_theta(blocks, len(design.y_lags), len(design.x)), design)).ravel()
    answer = lsmr(operator, residual, atol=1e-10, btol=1e-10, maxiter=min(size, 200))[0]
    directions: list[FactorBlock] = []
    for start, middle, end, width in layout:
        directions.append(
            FactorBlock(
                answer[start:middle].reshape(n, width),
                answer[middle:end].reshape(t, width),
            )
        )
    return directions


def _gauss_newton_refine(
    y: Array,
    design: Design,
    blocks: list[FactorBlock],
    coefficient_bound: float,
    interior_tolerance: float,
    current_objective: float,
) -> tuple[list[FactorBlock], float, bool]:
    directions = _gauss_newton_direction(y, design, blocks)
    step = 1.0
    for _ in range(12):
        candidate = [
            FactorBlock(
                block.loading + step * direction.loading,
                block.factor + step * direction.factor,
            )
            for block, direction in zip(blocks, directions, strict=True)
        ]
        candidate = rebalance_blocks(candidate)
        theta = _to_theta(candidate, len(design.y_lags), len(design.x))
        value = _objective(y, design, candidate)
        if (
            max_abs(theta) < coefficient_bound - interior_tolerance
            and value <= current_objective + 1e-11 * max(1.0, current_objective)
            and value < current_objective
        ):
            return candidate, value, True
        step *= 0.5
    return blocks, current_objective, False


def _fit_interior(
    y: Array,
    design: Design,
    widths: tuple[int, ...],
    *,
    seed: int | np.random.SeedSequence,
    max_sweeps: int,
    objective_rtol: float,
    stationarity_tol: float,
    coefficient_bound: float,
    interior_numerical_tolerance: float,
    lstsq_rcond: float,
    diagnostic_context: str | None,
) -> FitResult:
    started = time.perf_counter()
    rng = np.random.default_rng(seed)
    blocks = _initial_blocks(y.shape, widths, rng, None)
    regressors = design.regressors()
    initial_envelope = max_abs(
        _to_theta(blocks, len(design.y_lags), len(design.x))
    )
    history = [_objective(y, design, blocks)]
    stationarity_history: list[float] = []
    condition_history: list[float] = []
    safeguards = 0
    gn_steps = 0
    objective_stopped = False
    for _sweep in range(1, max_sweeps + 1):
        previous_blocks = _copy_blocks(blocks)
        previous_value = history[-1]
        conditions: list[float] = []
        for i in range(y.shape[0]):
            pieces = [
                regressors[j][i, :, None] * block.factor
                for j, block in enumerate(blocks)
                if block.loading.shape[1]
            ]
            if pieces:
                solution, condition = _stable_lstsq(
                    np.concatenate(pieces, axis=1), y[i], lstsq_rcond
                )
                conditions.append(condition)
                offset = 0
                for block in blocks:
                    width = block.loading.shape[1]
                    if width:
                        block.loading[i] = solution[offset : offset + width]
                        offset += width
        blocks = rebalance_blocks(blocks)
        for column in range(y.shape[1]):
            pieces = [
                regressors[j][:, column, None] * block.loading
                for j, block in enumerate(blocks)
                if block.factor.shape[1]
            ]
            if pieces:
                solution, condition = _stable_lstsq(
                    np.concatenate(pieces, axis=1), y[:, column], lstsq_rcond
                )
                conditions.append(condition)
                offset = 0
                for block in blocks:
                    width = block.factor.shape[1]
                    if width:
                        block.factor[column] = solution[offset : offset + width]
                        offset += width
        blocks = rebalance_blocks(blocks)
        value = _objective(y, design, blocks)
        if not np.isfinite(value) or value > previous_value + 1e-11 * max(1.0, previous_value):
            blocks = previous_blocks
            value = previous_value
            safeguards += 1
        theta = _to_theta(blocks, len(design.y_lags), len(design.x))
        stationarity = _stationarity(y, design, theta, widths)
        stationarity_history.append(stationarity)
        condition_history.append(max(conditions, default=1.0))
        history.append(value)
        relative = abs(previous_value - value) / max(abs(previous_value), 1e-14)
        if stationarity <= stationarity_tol and relative <= objective_rtol:
            objective_stopped = True
            break
        if relative <= objective_rtol and max_abs(theta) < coefficient_bound - interior_numerical_tolerance:
            refined, refined_value, accepted = _gauss_newton_refine(
                y,
                design,
                blocks,
                coefficient_bound,
                interior_numerical_tolerance,
                value,
            )
            if accepted:
                blocks = refined
                history[-1] = refined_value
                gn_steps += 1
            else:
                objective_stopped = True
                break
    theta = _to_theta(blocks, len(design.y_lags), len(design.x))
    stationarity = _stationarity(y, design, theta, widths)
    maximum = max_abs(theta)
    return FitResult(
        theta=theta,
        ranks=widths,
        objective=_objective(y, design, blocks),
        converged=bool(objective_stopped or stationarity <= stationarity_tol),
        iterations=len(history) - 1,
        objective_history=history,
        stationarity_residual=stationarity,
        max_envelope_ratio=maximum / coefficient_bound,
        factors=blocks,
        diagnostics={
            "solver_architecture": "gauge_aware_als_matrix_free_gauss_newton",
            "stationarity_pass": stationarity <= stationarity_tol,
            "stationarity_type": "product_tangent_projected_gradient",
            "stationarity_residual_history": stationarity_history,
            "condition_number_history": condition_history,
            "objective_safeguard_count": safeguards,
            "gauss_newton_steps": gn_steps,
            "bound_active": maximum >= coefficient_bound,
            "singular_values": singular_values(theta),
            "runtime_seconds": time.perf_counter() - started,
            "diagnostic_context": diagnostic_context,
            "initial_coefficient_envelope": initial_envelope,
            "final_coefficient_envelope": maximum,
        },
    )


def _fit_constrained(
    y: Array,
    design: Design,
    widths: tuple[int, ...],
    *,
    initial: Coefficients,
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
    initial_max = max_abs(initial)
    start = (
        scale(initial, coefficient_bound * (1.0 - 10.0 * constraint_tolerance) / initial_max)
        if initial_max >= coefficient_bound
        else initial
    )
    blocks = [
        balanced_factor_block(matrix, width)
        for matrix, width in zip(start.matrices(), widths, strict=True)
    ]
    regressors = design.regressors()
    history = [_objective(y, design, blocks)]
    kkt_history: list[float] = []
    condition_history: list[float] = []
    subproblem_iterations = 0
    solver_failed = False
    messages: list[str] = []
    loading_active: dict[int, tuple[int, ...]] = {}
    factor_active: dict[int, tuple[int, ...]] = {}
    kkt = float("inf")
    for _sweep in range(1, max_sweeps + 1):
        previous_blocks = _copy_blocks(blocks)
        previous_value = history[-1]
        sweep_conditions: list[float] = []
        for i in range(y.shape[0]):
            sub_design, current, constraints = _box_subproblem_matrices(
                blocks, regressors, i, update_loading=True
            )
            result = _active_set_qp(
                sub_design,
                y[i],
                constraints,
                current,
                coefficient_bound,
                max_iterations=constrained_subproblem_max_iterations,
                tolerance=constrained_subproblem_tolerance,
                constraint_tolerance=constraint_tolerance,
                warm_active_set=loading_active.get(i, ()),
            )
            subproblem_iterations += result.iterations
            sweep_conditions.append(result.condition_number)
            if not result.accepted:
                messages.append(f"loading[{i}]: {result.status}")
                solver_failed = True
                break
            loading_active[i] = result.active_set
            offset = 0
            for block in blocks:
                width = block.loading.shape[1]
                if width:
                    block.loading[i] = result.solution[offset : offset + width]
                    offset += width
        if solver_failed:
            blocks = previous_blocks
            break
        for block in blocks:
            _renormalize(block)
        for column in range(y.shape[1]):
            sub_design, current, constraints = _box_subproblem_matrices(
                blocks, regressors, column, update_loading=False
            )
            result = _active_set_qp(
                sub_design,
                y[:, column],
                constraints,
                current,
                coefficient_bound,
                max_iterations=constrained_subproblem_max_iterations,
                tolerance=constrained_subproblem_tolerance,
                constraint_tolerance=constraint_tolerance,
                warm_active_set=factor_active.get(column, ()),
            )
            subproblem_iterations += result.iterations
            sweep_conditions.append(result.condition_number)
            if not result.accepted:
                messages.append(f"factor[{column}]: {result.status}")
                solver_failed = True
                break
            factor_active[column] = result.active_set
            offset = 0
            for block in blocks:
                width = block.factor.shape[1]
                if width:
                    block.factor[column] = result.solution[offset : offset + width]
                    offset += width
        value = _objective(y, design, blocks)
        if solver_failed or value > previous_value + 1e-11 * max(1.0, previous_value):
            blocks = previous_blocks
            if not solver_failed:
                messages.append("full_sweep_objective_increase")
                solver_failed = True
            break
        history.append(value)
        kkt = _factor_kkt_residual(
            blocks,
            regressors,
            y,
            coefficient_bound,
            max(10.0 * constraint_tolerance, 1e-8),
        )
        kkt_history.append(kkt)
        condition_history.append(max(sweep_conditions, default=1.0))
        relative = abs(previous_value - value) / max(1.0, abs(previous_value))
        if relative <= objective_rtol and kkt <= constrained_kkt_tolerance:
            break
    theta = _to_theta(blocks, len(design.y_lags), len(design.x))
    maximum = max_abs(theta)
    violation = max(0.0, maximum - coefficient_bound)
    finite = bool(np.isfinite(_objective(y, design, blocks)) and all(np.all(np.isfinite(matrix)) for matrix in theta.matrices()))
    feasibility_pass = finite and violation <= constraint_tolerance
    if not np.isfinite(kkt):
        kkt = _factor_kkt_residual(
            blocks,
            regressors,
            y,
            coefficient_bound,
            max(10.0 * constraint_tolerance, 1e-8),
        )
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
    runtime = time.perf_counter() - started
    return FitResult(
        theta=theta,
        ranks=widths,
        objective=_objective(y, design, blocks),
        converged=converged,
        iterations=len(history) - 1,
        objective_history=history,
        stationarity_residual=kkt,
        max_envelope_ratio=maximum / coefficient_bound,
        factors=blocks,
        diagnostics={
            "solver_architecture": "reversibly_preconditioned_active_set_qp",
            "stationarity_pass": kkt_pass,
            "stationarity_type": "factor_space_box_KKT",
            "stationarity_residual_history": kkt_history,
            "condition_number_history": condition_history,
            "bound_active": maximum >= coefficient_bound - 10.0 * constraint_tolerance,
            "boundary_active": maximum >= coefficient_bound - 10.0 * constraint_tolerance,
            "singular_values": singular_values(theta),
            "runtime_seconds": runtime,
            "constrained_runtime": runtime,
            "diagnostic_context": diagnostic_context,
            "initial_coefficient_envelope": initial_max,
            "final_coefficient_envelope": maximum,
            "coefficient_envelope_history": [max_abs(start), maximum],
            "max_constraint_violation": violation,
            "constraint_tolerance": constraint_tolerance,
            "constrained_KKT_residual": kkt,
            "constrained_iterations": subproblem_iterations,
            "constrained_solver_status": status,
            "constrained_objective": _objective(y, design, blocks),
            "constrained_algorithm": "deterministic_active_set_linear_box_QP",
            "constrained_messages": messages,
        },
    )


def fit_cap_plus_one(
    y: Array,
    design: Design,
    widths: tuple[int, ...],
    *,
    seed: int | np.random.SeedSequence = 0,
    max_sweeps: int = 200,
    objective_rtol: float = 1e-8,
    stationarity_tol: float = 1e-6,
    coefficient_bound: float = 10.0,
    lstsq_rcond: float = 1e-10,
    interior_numerical_tolerance: float = 1e-8,
    constraint_tolerance: float = 1e-8,
    constrained_kkt_tolerance: float = 1e-4,
    constrained_subproblem_tolerance: float = 1e-10,
    constrained_subproblem_max_iterations: int = 200,
    diagnostic_context: str | None = None,
    **_unused: Any,
) -> FitResult:
    """Solve the frozen rank-at-most cap+1 coefficient problem."""

    if design.shape != y.shape or len(widths) != len(design.y_lags) + len(design.x) + 1:
        raise ValueError("incompatible data, design, and widths")
    interior = _fit_interior(
        y,
        design,
        widths,
        seed=seed,
        max_sweeps=max_sweeps,
        objective_rtol=objective_rtol,
        stationarity_tol=stationarity_tol,
        coefficient_bound=coefficient_bound,
        interior_numerical_tolerance=interior_numerical_tolerance,
        lstsq_rcond=lstsq_rcond,
        diagnostic_context=diagnostic_context,
    )
    unconstrained_max = max_abs(interior.theta)
    if unconstrained_max < coefficient_bound - interior_numerical_tolerance:
        interior.diagnostics.update(
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
        result = interior
    else:
        result = _fit_constrained(
            y,
            design,
            widths,
            initial=interior.theta,
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
                "unconstrained_objective": interior.objective,
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
