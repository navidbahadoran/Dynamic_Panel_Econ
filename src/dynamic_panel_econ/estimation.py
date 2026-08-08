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
    coefficient_bound: float = 9.0,
    lstsq_rcond: float = 1e-10,
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
        },
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
    coefficient_bound: float = 9.0,
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
