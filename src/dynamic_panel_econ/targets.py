"""Linear target directions and exact truth evaluation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import Coefficients, inner, zeros_like


@dataclass(frozen=True, slots=True)
class TargetSpec:
    name: str
    direction: Coefficients
    broad: bool


def paper_index(size: int, divisor: int) -> int:
    """Convert mathematical ``floor(size/divisor)`` (1-based) to Python indexing."""

    mathematical = max(1, size // divisor)
    return mathematical - 1


def _set_block(direction: Coefficients, block: str, weights: np.ndarray) -> None:
    if block == "A":
        direction.A[0][:] = weights
    elif block == "B":
        direction.B[0][:] = weights
    else:
        raise ValueError("target block must be A or B")


def target_direction(
    name: str, template: Coefficients, groups: np.ndarray | None = None
) -> TargetSpec:
    n, t = template.shape
    direction = zeros_like(template)
    weights = np.zeros((n, t), dtype=np.float64)
    i0, t0 = paper_index(n, 4), paper_index(t, 2)
    parts = name.split("_")
    block = parts[0]
    broad = False
    if name.endswith("_entry"):
        weights[i0, t0] = 1.0
    elif name.endswith("_fixed_time_mean"):
        weights[:, t0] = 1.0 / n
    elif name.endswith("_full_mean"):
        weights[:] = 1.0 / (n * t)
        broad = True
    elif any(token in name for token in ("G1", "G2")):
        if groups is None:
            raise ValueError("group target requires group assignments")
        g1, g2 = groups == 0, groups == 1
        fixed = "fixed_time" in name
        broad = not fixed
        if "G2_minus_G1" in name:
            if fixed:
                weights[g2, t0] = 1.0 / g2.sum()
                weights[g1, t0] = -1.0 / g1.sum()
            else:
                weights[g2, :] = 1.0 / (g2.sum() * t)
                weights[g1, :] = -1.0 / (g1.sum() * t)
        else:
            group = g1 if "G1" in name else g2
            if fixed:
                weights[group, t0] = 1.0 / group.sum()
            else:
                weights[group, :] = 1.0 / (group.sum() * t)
    else:
        raise ValueError(f"unknown target: {name}")
    _set_block(direction, block, weights)
    return TargetSpec(name, direction, broad)


def target_value(direction: Coefficients, theta: Coefficients) -> float:
    return inner(direction, theta)


def default_target_names(include_groups: bool) -> list[str]:
    names = [
        "A_entry",
        "B_entry",
        "A_fixed_time_mean",
        "B_fixed_time_mean",
        "A_full_mean",
        "B_full_mean",
    ]
    if include_groups:
        for block in ("A", "B"):
            for target in ("G1", "G2", "G2_minus_G1"):
                names.append(f"{block}_{target}_fixed_time")
                names.append(f"{block}_{target}_time_average")
    return names


def restrict_direction(
    direction: Coefficients, rows: np.ndarray, cols: np.ndarray
) -> Coefficients:
    ix = np.ix_(rows, cols)
    return Coefficients(
        [matrix[ix] for matrix in direction.A],
        [matrix[ix] for matrix in direction.B],
        direction.H[ix],
    )


def embedded_restriction(
    direction: Coefficients, rows: np.ndarray, cols: np.ndarray
) -> Coefficients:
    output = zeros_like(direction)
    ix = np.ix_(rows, cols)
    for target, source in zip(output.matrices(), direction.matrices(), strict=True):
        target[ix] = source[ix]
    return output
