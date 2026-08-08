"""Coefficient-collection algebra and the dynamic-panel design operator."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


@dataclass(slots=True)
class Coefficients:
    """Ordered low-rank coefficient matrices for lags, covariates, and intercept."""

    A: list[Array]
    B: list[Array]
    H: Array

    @property
    def shape(self) -> tuple[int, int]:
        return self.H.shape

    def matrices(self) -> list[Array]:
        return [*self.A, *self.B, self.H]

    def copy(self) -> Coefficients:
        return Coefficients([x.copy() for x in self.A], [x.copy() for x in self.B], self.H.copy())


@dataclass(slots=True)
class Design:
    """Observed regressors aligned to response times, all with shape ``(N,T)``."""

    y_lags: list[Array]
    x: list[Array]

    @property
    def shape(self) -> tuple[int, int]:
        arrays = [*self.y_lags, *self.x]
        if not arrays:
            raise ValueError("a design must contain at least one regressor")
        return arrays[0].shape

    def regressors(self, include_intercept: bool = True) -> list[Array]:
        out = [*self.y_lags, *self.x]
        if include_intercept:
            out.append(np.ones(self.shape, dtype=np.float64))
        return out


def validate_collection(theta: Coefficients, design: Design) -> None:
    if len(theta.A) != len(design.y_lags) or len(theta.B) != len(design.x):
        raise ValueError("coefficient and design block counts differ")
    if any(m.shape != design.shape for m in theta.matrices()):
        raise ValueError("all coefficient and design matrices must share shape (N,T)")


def fitted_values(theta: Coefficients, design: Design) -> Array:
    """Apply the exact fitted-value map to a coefficient collection."""

    validate_collection(theta, design)
    fitted = theta.H.copy()
    for matrix, regressor in zip(theta.A, design.y_lags, strict=True):
        fitted += matrix * regressor
    for matrix, regressor in zip(theta.B, design.x, strict=True):
        fitted += matrix * regressor
    return fitted


def adjoint(residual: Array, design: Design) -> Coefficients:
    """Apply the Frobenius adjoint of :func:`fitted_values`."""

    if residual.shape != design.shape:
        raise ValueError("residual and design shapes differ")
    return Coefficients(
        [residual * z for z in design.y_lags],
        [residual * z for z in design.x],
        residual.copy(),
    )


def zeros_like(theta: Coefficients) -> Coefficients:
    return Coefficients(
        [np.zeros_like(x) for x in theta.A],
        [np.zeros_like(x) for x in theta.B],
        np.zeros_like(theta.H),
    )


def zeros_for_design(design: Design) -> Coefficients:
    n, t = design.shape
    return Coefficients(
        [np.zeros((n, t), dtype=np.float64) for _ in design.y_lags],
        [np.zeros((n, t), dtype=np.float64) for _ in design.x],
        np.zeros((n, t), dtype=np.float64),
    )


def inner(left: Coefficients, right: Coefficients) -> float:
    if len(left.A) != len(right.A) or len(left.B) != len(right.B):
        raise ValueError("coefficient block counts differ")
    return float(sum(np.vdot(x, y) for x, y in zip(left.matrices(), right.matrices(), strict=True)))


def add(left: Coefficients, right: Coefficients, alpha: float = 1.0) -> Coefficients:
    mats = [x + alpha * y for x, y in zip(left.matrices(), right.matrices(), strict=True)]
    p, k = len(left.A), len(left.B)
    return Coefficients(mats[:p], mats[p : p + k], mats[-1])


def scale(theta: Coefficients, value: float) -> Coefficients:
    mats = [value * x for x in theta.matrices()]
    p, k = len(theta.A), len(theta.B)
    return Coefficients(mats[:p], mats[p : p + k], mats[-1])


def frobenius_norm(theta: Coefficients) -> float:
    return float(np.sqrt(sum(np.vdot(x, x) for x in theta.matrices())))


def max_abs(theta: Coefficients) -> float:
    return float(max(np.max(np.abs(x), initial=0.0) for x in theta.matrices()))


def from_matrices(matrices: Iterable[Array], p: int, k: int) -> Coefficients:
    mats = list(matrices)
    if len(mats) != p + k + 1:
        raise ValueError("wrong number of coefficient matrices")
    return Coefficients(mats[:p], mats[p : p + k], mats[-1])


def subset_design(design: Design, rows: NDArray[np.int_], cols: NDArray[np.int_]) -> Design:
    ix = np.ix_(rows, cols)
    return Design([z[ix] for z in design.y_lags], [z[ix] for z in design.x])


def subset_coefficients(
    theta: Coefficients, rows: NDArray[np.int_], cols: NDArray[np.int_]
) -> Coefficients:
    ix = np.ix_(rows, cols)
    return Coefficients([m[ix] for m in theta.A], [m[ix] for m in theta.B], theta.H[ix])
