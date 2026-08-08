"""One-dimensional positive-semidefinite Bartlett spatial HAC."""

from __future__ import annotations

import numpy as np


def spatial_cutoff(n: int, t: int, c_sp: float = 1.0) -> int:
    return int(np.ceil(c_sp * np.log(n * t)))


def bartlett_quadratic(scores: np.ndarray, cutoff: int) -> float:
    """Compute the Bartlett quadratic form in O(T N h), without dense weights."""

    n, _ = scores.shape
    h = min(max(int(cutoff), 0), n - 1)
    value = float(np.sum(scores * scores))
    for distance in range(1, h + 1):
        weight = 1.0 - distance / (h + 1.0)
        value += 2.0 * weight * float(np.sum(scores[:-distance] * scores[distance:]))
    return value


def bartlett_matrix(n: int, cutoff: int) -> np.ndarray:
    h = min(max(int(cutoff), 0), n - 1)
    index = np.arange(n)
    distances = np.abs(index[:, None] - index[None, :])
    return np.maximum(1.0 - distances / (h + 1.0), 0.0)


def bartlett_quadratic_dense(scores: np.ndarray, cutoff: int) -> float:
    weights = bartlett_matrix(scores.shape[0], cutoff)
    return float(sum(scores[:, column] @ weights @ scores[:, column] for column in range(scores.shape[1])))
