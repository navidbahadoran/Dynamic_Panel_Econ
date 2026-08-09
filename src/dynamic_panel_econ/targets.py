"""Linear target directions and exact truth evaluation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import Coefficients, inner, zeros_like
from .lowrank import numerical_rank, tangent_project


@dataclass(frozen=True, slots=True)
class TargetSpec:
    name: str
    direction: Coefficients
    broad: bool
    theorem_validation: bool = True
    applicability: str = "theorem_covered"


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
    name: str,
    template: Coefficients,
    groups: np.ndarray | None = None,
    *,
    dgp: int | None = None,
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
    theorem_validation = True
    applicability = "theorem_covered"
    if "G2_minus_G1_fixed_time" in name and dgp is not None and dgp < 4:
        theorem_validation = False
        applicability = "weak_target_stress_outside_assumption9"
    return TargetSpec(name, direction, broad, theorem_validation, applicability)


def target_value(direction: Coefficients, theta: Coefficients) -> float:
    return inner(direction, theta)


def target_regularity_diagnostics(
    spec: TargetSpec,
    theta0: Coefficients,
) -> dict[str, float | str | bool | None]:
    """Compute true tangent projection and entry leverage diagnostics."""

    ranks = tuple(numerical_rank(matrix) for matrix in theta0.matrices())
    direction_norm = float(
        np.sqrt(sum(np.vdot(matrix, matrix) for matrix in spec.direction.matrices()))
    )
    projected = tangent_project(spec.direction, theta0, ranks)
    projected_norm = float(
        np.sqrt(sum(np.vdot(matrix, matrix) for matrix in projected.matrices()))
    )
    diagnostics: dict[str, float | str | bool | None] = {
        "target_applicability": spec.applicability,
        "headline_theorem_target": spec.theorem_validation,
        "true_target_direction_norm": direction_norm,
        "true_target_tangent_norm": projected_norm,
        "true_target_projection_ratio": (
            projected_norm / direction_norm if direction_norm > 0.0 else float("nan")
        ),
        "true_entry_unit_leverage_scaled": None,
        "true_entry_time_leverage_scaled": None,
    }
    if spec.name.endswith("_entry"):
        block = theta0.A[0] if spec.name.startswith("A_") else theta0.B[0]
        rank = ranks[0] if spec.name.startswith("A_") else ranks[1]
        i0, t0 = paper_index(block.shape[0], 4), paper_index(block.shape[1], 2)
        if rank > 0:
            u, _, vt = np.linalg.svd(block, full_matrices=False)
            diagnostics["true_entry_unit_leverage_scaled"] = float(
                block.shape[0] * np.sum(u[i0, :rank] ** 2)
            )
            diagnostics["true_entry_time_leverage_scaled"] = float(
                block.shape[1] * np.sum(vt[:rank, t0] ** 2)
            )
        else:
            diagnostics["true_entry_unit_leverage_scaled"] = 0.0
            diagnostics["true_entry_time_leverage_scaled"] = 0.0
    return diagnostics


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
