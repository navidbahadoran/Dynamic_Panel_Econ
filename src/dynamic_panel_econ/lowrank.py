"""Low-rank matrix primitives shared by estimation and inference."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .core import Coefficients, from_matrices

Array = NDArray[np.float64]


def compact_svd(matrix: Array, rank: int | None = None) -> tuple[Array, Array, Array]:
    u, singular, vt = np.linalg.svd(matrix, full_matrices=False)
    if rank is not None:
        u, singular, vt = u[:, :rank], singular[:rank], vt[:rank]
    return u, singular, vt


def truncated_matrix(matrix: Array, rank: int) -> Array:
    if rank == 0:
        return np.zeros_like(matrix)
    u, singular, vt = compact_svd(matrix, rank)
    return (u * singular) @ vt


def threshold_rank(matrix: Array, threshold: float, cap: int) -> int:
    singular = np.linalg.svd(matrix, compute_uv=False)
    return int(min(cap, np.count_nonzero(singular > threshold)))


def tangent_project_matrix(matrix: Array, fitted: Array, rank: int) -> Array:
    """Project onto the exact-rank tangent space at ``fitted``."""

    if rank == 0:
        return np.zeros_like(matrix)
    u, _, vt = compact_svd(fitted, rank)
    v = vt.T
    left = u @ (u.T @ matrix)
    right = (matrix @ v) @ v.T
    overlap = u @ ((u.T @ matrix) @ v) @ v.T
    return left + right - overlap


def tangent_project(theta: Coefficients, fitted: Coefficients, ranks: tuple[int, ...]) -> Coefficients:
    matrices = [
        tangent_project_matrix(z, m, rank)
        for z, m, rank in zip(theta.matrices(), fitted.matrices(), ranks, strict=True)
    ]
    return from_matrices(matrices, len(theta.A), len(theta.B))


def numerical_rank(matrix: Array, tolerance: float | None = None) -> int:
    singular = np.linalg.svd(matrix, compute_uv=False)
    if not singular.size:
        return 0
    threshold = tolerance if tolerance is not None else max(matrix.shape) * np.finfo(float).eps * singular[0]
    return int(np.count_nonzero(singular > threshold))


def singular_values(theta: Coefficients) -> list[list[float]]:
    return [np.linalg.svd(matrix, compute_uv=False).tolist() for matrix in theta.matrices()]
