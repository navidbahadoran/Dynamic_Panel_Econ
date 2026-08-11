"""Revision-10 ridge ratios and preserved legacy Revision-9 IC selection."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .core import (
    Coefficients,
    Design,
    adjoint,
    fitted_values,
    from_matrices,
    max_abs,
    scale,
    zeros_like,
)
from .estimation import FitResult, NuclearFit, adapt_initial, fit_fixed_rank, nuclear_path
from .lowrank import numerical_rank, tangent_project, threshold_rank, truncated_matrix

RankVector = tuple[int, ...]
BEST_BASIN_PERTURBATION_MAGNITUDES = (1e-4, 3e-4, 1e-3)


class RankSelectionFailure(RuntimeError):
    """Raised when no valid post-refit can determine the IC minimum."""


class RankPilotFailure(RankSelectionFailure):
    """Raised when the imposed rank-cap pilot fails a paper diagnostic."""

    def __init__(self, message: str, diagnostics: dict[str, Any] | None = None) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics or {}


class FinalPostRefitFailure(RankSelectionFailure):
    """Raised when the selected-rank literal post-refit is numerically unresolved."""

    def __init__(self, message: str, diagnostics: dict[str, Any] | None = None) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics or {}


@dataclass(slots=True)
class Revision10RankSelectionResult:
    """Frozen Revision-10 spectral-pilot selection and final post-refit."""

    selected_ranks: RankVector
    pilot_fit: FitResult
    final_fit: FitResult
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class CandidateFit:
    ranks: RankVector
    fit: FitResult
    ic: float
    dimension: int
    sources: list[str]
    third_start_used: bool
    valid: bool
    invalid_reasons: list[str]
    start_diagnostics: dict[str, Any]


@dataclass(slots=True)
class RankSelectionResult:
    selected: CandidateFit
    candidates: dict[RankVector, CandidateFit]
    nuclear_fits: list[NuclearFit]
    cap_fit: FitResult
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class CapPilotRoute:
    """One deterministic, distinct-rank nuclear-space cap-pilot route."""

    path_index: int | None
    rank: RankVector
    theta: Coefficients
    penalty: float | None
    nuclear_objective: float | None
    first_path_index: int | None
    best_path_index: int | None
    source: str


def model_dimension(ranks: RankVector, n: int, t: int) -> int:
    return int(sum(rank * (n + t - rank) for rank in ranks))


def revision8_kappa(
    n: int,
    t: int,
    *,
    eta_for_penalty: float = 4.0,
    spatial_dimension: int = 1,
    multiplier: float = 1.0,
) -> float:
    """Revision-8 IC penalty: ``b_NT^2 log(NT)^(d_s+3)``."""

    nt = n * t
    log_nt = np.log(nt)
    b_nt = nt ** (1.0 / (8.0 + eta_for_penalty)) * log_nt
    return float(multiplier * b_nt**2 * log_nt ** (spatial_dimension + 3))


def information_criterion(
    objective: float,
    ranks: RankVector,
    n: int,
    t: int,
    eta_for_penalty: float = 4.0,
    spatial_dimension: int = 1,
    multiplier: float = 1.0,
) -> float:
    qhat = max(2.0 * objective, np.finfo(float).tiny)
    penalty = revision8_kappa(
        n,
        t,
        eta_for_penalty=eta_for_penalty,
        spatial_dimension=spatial_dimension,
        multiplier=multiplier,
    )
    return float(np.log(qhat) + penalty * model_dimension(ranks, n, t) / (n * t))


def one_coordinate_neighbors(ranks: RankVector, caps: RankVector) -> set[RankVector]:
    neighbors: set[RankVector] = set()
    for index, (rank, cap) in enumerate(zip(ranks, caps, strict=True)):
        if rank > 0:
            candidate = list(ranks)
            candidate[index] -= 1
            neighbors.add(tuple(candidate))
        if rank < cap:
            candidate = list(ranks)
            candidate[index] += 1
            neighbors.add(tuple(candidate))
    return neighbors


def _rank_vector(theta: Coefficients, threshold: float, caps: RankVector) -> RankVector:
    return tuple(
        threshold_rank(matrix, threshold, cap)
        for matrix, cap in zip(theta.matrices(), caps, strict=True)
    )


def _numerical_rank_vector(theta: Coefficients) -> RankVector:
    return tuple(numerical_rank(matrix) for matrix in theta.matrices())


def fit_invalid_reasons(
    fit: FitResult,
    required_ranks: RankVector,
    stationarity_tolerance: float,
    *,
    require_exact_numerical_rank: bool = True,
) -> list[str]:
    reasons = []
    if not fit.converged:
        reasons.append(str(fit.diagnostics.get("constrained_solver_status", "not_converged")))
    if fit.diagnostics.get("constrained_fallback_used"):
        if not bool(fit.diagnostics.get("stationarity_pass", False)):
            reasons.append("constrained_optimality_failure")
        if float(fit.diagnostics.get("max_constraint_violation", np.inf)) > float(
            fit.diagnostics.get("constraint_tolerance", 0.0)
        ):
            reasons.append("constrained_feasibility_failure")
    elif (
        not np.isfinite(fit.stationarity_residual)
        or fit.stationarity_residual > stationarity_tolerance
    ):
        reasons.append("stationarity_high")
    if not np.isfinite(fit.max_envelope_ratio) or fit.max_envelope_ratio > 1.0 + 1e-8:
        reasons.append("constrained_feasibility_failure")
    actual = _numerical_rank_vector(fit.theta)
    if require_exact_numerical_rank and actual != required_ranks:
        reasons.append(f"numerical_rank_support:{actual}")
    if not np.isfinite(fit.objective):
        reasons.append("nonfinite_objective")
    return reasons


def revision10_block_names(design: Design) -> tuple[str, ...]:
    """Ordered coefficient-block names used by Revision-10 Section 4.5."""

    return tuple(
        [*(f"A{index}" for index in range(1, len(design.y_lags) + 1)),
         *(f"B{index}" for index in range(1, len(design.x) + 1)), "H"]
    )


def revision10_scale_weights(design: Design) -> tuple[float, ...]:
    """Return the uncentered full-sample RMS weights in Revision-10 Section 4.5."""

    weights = [
        *(float(np.sqrt(np.mean(np.square(regressor)))) for regressor in design.y_lags),
        *(float(np.sqrt(np.mean(np.square(regressor)))) for regressor in design.x),
        1.0,
    ]
    if not design.y_lags:
        raise ValueError("Revision-10 rank selection requires reference block A1")
    if not all(np.isfinite(value) and value >= 0.0 for value in weights):
        raise RankPilotFailure("rank_selection_numerically_unresolved: nonfinite scale weight")
    reference = weights[0]
    if not np.isfinite(reference) or reference <= 0.0:
        raise RankPilotFailure(
            "rank_selection_numerically_unresolved: unusable reference weight w_A1",
            {"scale_weights": weights, "reference_weight_w_A1": reference},
        )
    return tuple(weights)


def revision10_normalized_spectrum(
    matrix: np.ndarray,
    *,
    block_weight: float,
    reference_weight: float,
    n: int,
    t: int,
    count: int,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Compute the frozen normalized squared pilot spectrum through cap+1."""

    if count < 1 or count > min(matrix.shape):
        raise ValueError("requested spectrum count must be between one and min(N,T)")
    if matrix.shape != (n, t):
        raise ValueError("pilot matrix and panel dimensions differ")
    if not np.isfinite(reference_weight) or reference_weight <= 0.0:
        raise RankPilotFailure(
            "rank_selection_numerically_unresolved: unusable reference weight w_A1"
        )
    singular = np.linalg.svd(matrix, compute_uv=False)[:count]
    scale_factor = (block_weight / reference_weight) ** 2 / (n * t)
    normalized = scale_factor * np.square(singular)
    if not np.all(np.isfinite(normalized)):
        raise RankPilotFailure(
            "rank_selection_numerically_unresolved: nonfinite normalized pilot spectrum"
        )
    return tuple(float(value) for value in singular), tuple(
        float(value) for value in normalized
    )


def revision10_ridge(n: int, t: int) -> float:
    """Return exactly ``a_NT=1/log(NT)`` from Revision-10 Section 4.5."""

    if n * t <= 1:
        raise ValueError("Revision-10 ridge requires N*T > 1")
    return float(1.0 / np.log(n * t))


def revision10_ridge_ratios(
    normalized_spectrum: tuple[float, ...] | list[float] | np.ndarray,
    *,
    reporting_cap: int,
    n: int,
    t: int,
) -> tuple[float, ...]:
    """Form the rank-zero anchor and positive-rank ratios, including cap+1."""

    spectrum = np.asarray(normalized_spectrum, dtype=float)
    if reporting_cap < 0 or spectrum.size < reporting_cap + 1:
        raise ValueError("normalized spectrum must include the genuine cap+1 value")
    if not np.all(np.isfinite(spectrum[: reporting_cap + 1])):
        raise RankPilotFailure(
            "rank_selection_numerically_unresolved: nonfinite ridge-ratio input"
        )
    ridge = revision10_ridge(n, t)
    ratios = [(spectrum[0] + ridge) / (1.0 + ridge)]
    ratios.extend(
        (spectrum[index] + ridge) / (spectrum[index - 1] + ridge)
        for index in range(1, reporting_cap + 1)
    )
    return tuple(float(value) for value in ratios)


def revision10_select_block_rank(ratios: tuple[float, ...] | list[float]) -> int:
    """Select the first exact argmin, so exact ties go to the smaller rank."""

    values = np.asarray(ratios, dtype=float)
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise RankPilotFailure(
            "rank_selection_numerically_unresolved: invalid ridge-ratio vector"
        )
    return int(np.argmin(values))


def revision10_assemble_rank_vector(
    block_ratios: list[tuple[float, ...]] | tuple[tuple[float, ...], ...],
) -> RankVector:
    """Assemble separately selected A, B, and H block ranks."""

    return tuple(revision10_select_block_rank(ratios) for ratios in block_ratios)


def build_candidates(
    preliminary: list[NuclearFit],
    cap_fit: FitResult,
    caps: RankVector,
    threshold: float,
) -> tuple[set[RankVector], dict[RankVector, list[str]], list[RankVector]]:
    sources: dict[RankVector, list[str]] = {}
    path_ranks = []
    for index, fit in enumerate(preliminary):
        rank = _rank_vector(fit.theta, threshold, caps)
        path_ranks.append(rank)
        sources.setdefault(rank, []).append(f"nuclear_path_{index}")
    cap_rank = _rank_vector(cap_fit.theta, threshold, caps)
    sources.setdefault(cap_rank, []).append("rank_cap_pilot")
    base = set(sources)
    candidates = set(base)
    for rank in base:
        for neighbor in one_coordinate_neighbors(rank, caps):
            candidates.add(neighbor)
            sources.setdefault(neighbor, []).append(f"neighbor_of_{rank}")
    return candidates, sources, path_ranks


def _closest_preliminary(
    ranks: RankVector, preliminary: list[NuclearFit], caps: RankVector, threshold: float
) -> Coefficients:
    scored = []
    for index, fit in enumerate(preliminary):
        proposal = _rank_vector(fit.theta, threshold, caps)
        distance = sum(abs(a - b) for a, b in zip(ranks, proposal, strict=True))
        scored.append((distance, index, fit.theta))
    return min(scored, key=lambda item: (item[0], item[1]))[2]


def _theta_fingerprint(theta: Coefficients) -> str:
    digest = hashlib.sha256()
    for matrix in theta.matrices():
        digest.update(np.ascontiguousarray(matrix).view(np.uint8))
    return digest.hexdigest()


def _rescale_cap_start(
    theta: Coefficients,
    coefficient_bound: float,
    start_envelope_fraction: float,
) -> tuple[Coefficients, dict[str, float | bool]]:
    """Common-scale one numerical start into a strict interior envelope."""

    original = max_abs(theta)
    applied = (
        min(1.0, start_envelope_fraction * coefficient_bound / original)
        if original > 0.0
        else 1.0
    )
    rescaled = scale(theta, applied)
    return rescaled, {
        "original_max_abs_coefficient": original,
        "applied_common_scale": applied,
        "rescaled_max_abs_coefficient": max_abs(rescaled),
        "common_rescaling_used": applied < 1.0,
    }


def _objective_at(y: np.ndarray, design: Design, theta: Coefficients) -> float:
    residual = y - fitted_values(theta, design)
    return float(np.vdot(residual, residual) / (2.0 * y.size))


def _rank_preserving_perturbation(
    theta: Coefficients,
    ranks: RankVector,
    magnitude: float,
) -> Coefficients:
    """Apply a deterministic canonical-SVD perturbation without clipping entries."""

    perturbed: list[np.ndarray] = []
    for block_index, (matrix, rank) in enumerate(
        zip(theta.matrices(), ranks, strict=True)
    ):
        if rank == 0:
            perturbed.append(np.zeros_like(matrix))
            continue
        u, singular, vt = np.linalg.svd(matrix, full_matrices=False)
        u = u[:, :rank]
        v = vt[:rank].T
        singular = singular[:rank]
        row_grid = np.arange(1, matrix.shape[0] + 1, dtype=float)[:, None]
        col_grid = np.arange(1, matrix.shape[1] + 1, dtype=float)[:, None]
        component_grid = np.arange(1, rank + 1, dtype=float)[None, :]
        du = np.sin((block_index + 1) * row_grid * component_grid)
        dv = np.cos((block_index + 2) * col_grid * component_grid)
        du -= u @ (u.T @ du)
        dv -= v @ (v.T @ dv)
        u_new, _ = np.linalg.qr(u + magnitude * du, mode="reduced")
        v_new, _ = np.linalg.qr(v + magnitude * dv, mode="reduced")
        signs = np.where(np.arange(rank) % 2 == 0, 1.0, -1.0)
        scale = max(float(singular[0]), 1.0)
        singular_new = singular + magnitude * scale * signs
        floor = max(np.finfo(float).eps * scale, magnitude**2 * scale)
        singular_new = np.maximum(singular_new, floor)
        perturbed.append((u_new * singular_new) @ v_new.T)
    return from_matrices(perturbed, len(theta.A), len(theta.B))


def _confirmation_record(
    fit: FitResult,
    start: Coefficients,
    base: Coefficients,
    start_id: str,
    magnitude: float,
    best_objective: float,
    ranks: RankVector,
    stationarity_tolerance: float,
    objective_tolerance: float,
    y: np.ndarray,
    design: Design,
) -> dict[str, Any]:
    reasons = fit_invalid_reasons(fit, ranks, stationarity_tolerance)
    gap = abs(fit.objective - best_objective) / max(1.0, abs(best_objective))
    return {
        "route_type": "basin_confirmation",
        "confirmation_start_id": start_id,
        "perturbation_magnitude": magnitude,
        "perturbation_norm": float(
            np.sqrt(
                sum(
                    np.linalg.norm(after - before) ** 2
                    for before, after in zip(
                        start.matrices(), base.matrices(), strict=True
                    )
                )
            )
        ),
        "starting_objective": _objective_at(y, design, start),
        "final_objective": fit.objective,
        "convergence": fit.converged,
        "stationarity": fit.stationarity_residual,
        "coefficient_envelope": max_abs(fit.theta),
        "numerical_rank": _numerical_rank_vector(fit.theta),
        "runtime": fit.diagnostics.get("runtime_seconds"),
        "objective_gap_to_best": gap,
        "valid": not reasons,
        "invalid_reasons": reasons,
        "confirmation_pass": not reasons and gap <= objective_tolerance,
    }


def _confirm_best_basin(
    y: np.ndarray,
    design: Design,
    best_fit: FitResult,
    ranks: RankVector,
    *,
    fit_options: dict[str, Any],
    seed: int | np.random.SeedSequence,
    stationarity_tolerance: float,
    objective_tolerance: float,
    start_envelope_fraction: float,
    diagnostic_context: str,
) -> tuple[list[FitResult], list[dict[str, Any]], bool]:
    """Try to reproduce one best basin from independent deterministic perturbations."""

    coefficient_bound = float(fit_options.get("coefficient_bound", 10.0))
    base, _ = _rescale_cap_start(
        adapt_initial(best_fit.theta, ranks), coefficient_bound, start_envelope_fraction
    )
    fits: list[FitResult] = []
    records: list[dict[str, Any]] = []
    for index, magnitude in enumerate(BEST_BASIN_PERTURBATION_MAGNITUDES, start=1):
        start = _rank_preserving_perturbation(base, ranks, magnitude)
        if max_abs(start) >= coefficient_bound or _numerical_rank_vector(start) != ranks:
            records.append(
                {
                    "route_type": "basin_confirmation",
                    "confirmation_start_id": f"confirmation_{index}",
                    "perturbation_magnitude": magnitude,
                    "perturbation_norm": float("nan"),
                    "starting_objective": _objective_at(y, design, start),
                    "final_objective": float("nan"),
                    "convergence": False,
                    "stationarity": float("nan"),
                    "coefficient_envelope": max_abs(start),
                    "numerical_rank": _numerical_rank_vector(start),
                    "runtime": 0.0,
                    "objective_gap_to_best": float("nan"),
                    "valid": False,
                    "invalid_reasons": ["confirmation_start_not_strictly_interior"],
                    "confirmation_pass": False,
                }
            )
            continue
        fit = fit_fixed_rank(
            y,
            design,
            ranks,
            initial=start,
            seed=seed,
            diagnostic_context=f"{diagnostic_context}:confirmation_{index}",
            **fit_options,
        )
        fits.append(fit)
        records.append(
            _confirmation_record(
                fit,
                start,
                base,
                f"confirmation_{index}",
                magnitude,
                best_fit.objective,
                ranks,
                stationarity_tolerance,
                objective_tolerance,
                y,
                design,
            )
        )
    return fits, records, sum(record["confirmation_pass"] for record in records) >= 2


def fit_fixed_rank_multistart(
    y: np.ndarray,
    design: Design,
    ranks: RankVector,
    *,
    seed: int | np.random.SeedSequence,
    fit_options: dict[str, Any],
    stationarity_tolerance: float,
    start_objective_stability_tol: float,
    start_envelope_fraction: float = 0.8,
) -> tuple[FitResult, dict[str, Any]]:
    """Verify a supplied-rank solution with three deterministic starts and confirmation."""

    rng = np.random.default_rng(seed)
    additional_seeds = [int(value) for value in rng.integers(0, 2**32 - 1, size=2)]
    start_seeds: list[int | np.random.SeedSequence] = [seed, *additional_seeds]
    original_fits = [
        fit_fixed_rank(
            y,
            design,
            ranks,
            seed=start_seed,
            diagnostic_context=f"full_fixed_rank:original_start_{index}",
            **fit_options,
        )
        for index, start_seed in enumerate(start_seeds, start=1)
    ]
    original_records = []
    for index, (start_seed, fit) in enumerate(
        zip(start_seeds, original_fits, strict=True), start=1
    ):
        reasons = fit_invalid_reasons(fit, ranks, stationarity_tolerance)
        original_records.append(
            {
                "route_type": "original",
                "start_id": f"original_start_{index}",
                "deterministic_seed": (
                    f"baseline:{start_seed.entropy}:{start_seed.spawn_key}"
                    if isinstance(start_seed, np.random.SeedSequence)
                    else str(start_seed)
                ),
                "final_objective": fit.objective,
                "convergence": fit.converged,
                "stationarity": fit.stationarity_residual,
                "coefficient_envelope": max_abs(fit.theta),
                "numerical_rank": _numerical_rank_vector(fit.theta),
                "runtime": fit.diagnostics.get("runtime_seconds"),
                "valid": not reasons,
                "invalid_reasons": reasons,
            }
        )
    valid_original = sorted(
        [fit for fit in original_fits if not fit_invalid_reasons(fit, ranks, stationarity_tolerance)],
        key=lambda fit: fit.objective,
    )
    credible = sorted(
        original_fits,
        key=lambda fit: (
            len(
                [
                    reason
                    for reason in fit_invalid_reasons(fit, ranks, stationarity_tolerance)
                    if reason != "constrained_feasibility_failure"
                ]
            ),
            fit.objective,
        ),
    )
    numerical_best_objective = min(fit.objective for fit in original_fits)
    for record in original_records:
        record["objective_gap_to_best"] = abs(
            float(record["final_objective"]) - numerical_best_objective
        ) / max(1.0, abs(numerical_best_objective))
    best = valid_original[0] if valid_original else credible[0]
    matching_original = [
        fit
        for fit in valid_original
        if abs(fit.objective - best.objective) / max(1.0, abs(best.objective))
        <= start_objective_stability_tol
    ]
    confirmation_fits: list[FitResult] = []
    confirmation_records: list[dict[str, Any]] = []
    confirmation_pass = False
    acceptance_basis = "original_route_stability"
    stable = len(matching_original) >= 2
    if not stable and valid_original:
        confirmation_fits, confirmation_records, confirmation_pass = _confirm_best_basin(
            y,
            design,
            best,
            ranks,
            fit_options=fit_options,
            seed=start_seeds[0],
            stationarity_tolerance=stationarity_tolerance,
            objective_tolerance=start_objective_stability_tol,
            start_envelope_fraction=start_envelope_fraction,
            diagnostic_context="full_fixed_rank",
        )
        stable = confirmation_pass
        acceptance_basis = "confirmed_best_basin" if stable else "failure"
    elif not stable:
        acceptance_basis = "failure"
    valid_confirmation = [
        fit
        for fit in confirmation_fits
        if not fit_invalid_reasons(fit, ranks, stationarity_tolerance)
    ]
    eligible = [*valid_original, *valid_confirmation]
    chosen = min(eligible, key=lambda fit: fit.objective) if eligible else best
    ordered_objectives = sorted(fit.objective for fit in valid_original)
    diagnostics = {
        "algorithm": "fixed_rank_three_deterministic_starts_with_best_basin_confirmation",
        "requested_rank": ranks,
        "objective_stability_tolerance": start_objective_stability_tol,
        "original_start_count": 3,
        "original_start_records": original_records,
        "original_best_objective": ordered_objectives[0] if ordered_objectives else best.objective,
        "original_second_best_objective": (
            ordered_objectives[1] if len(ordered_objectives) >= 2 else float("nan")
        ),
        "original_stability_gap": (
            abs(ordered_objectives[1] - ordered_objectives[0])
            / max(1.0, abs(ordered_objectives[0]))
            if len(ordered_objectives) >= 2
            else float("nan")
        ),
        "confirmation_start_records": confirmation_records,
        "confirmation_best_objective": min(
            (record["final_objective"] for record in confirmation_records if record["valid"]),
            default=float("nan"),
        ),
        "number_confirmation_valid": sum(record["valid"] for record in confirmation_records),
        "number_confirmation_matching_best": sum(
            record["confirmation_pass"] for record in confirmation_records
        ),
        "objective_stability_pass": stable,
        "final_acceptance_basis": acceptance_basis,
    }
    return chosen, diagnostics


def _synthetic_cap_start(
    template: Coefficients,
    ranks: RankVector,
    threshold: float,
) -> Coefficients:
    """Create a deterministic, very-low-rank singular-space start."""

    matrices = []
    for matrix, rank in zip(template.matrices(), ranks, strict=True):
        proposal = np.zeros_like(matrix)
        for component in range(rank):
            row = np.cos((component + 1) * np.arange(matrix.shape[0], dtype=float))
            col = np.cos((component + 1) * np.arange(matrix.shape[1], dtype=float))
            row /= np.linalg.norm(row)
            col /= np.linalg.norm(col)
            proposal += (1.25 + 0.05 * component) * threshold * np.outer(row, col)
        matrices.append(proposal)
    return from_matrices(matrices, len(template.A), len(template.B))


def _cap_pilot_routes(
    preliminary: list[NuclearFit],
    threshold: float,
    caps: RankVector,
    *,
    max_routes: int,
    coefficient_bound: float,
) -> tuple[list[CapPilotRoute], list[dict[str, Any]]]:
    """Select bounded, distinct-rank routes from the complete nuclear path."""

    occurrences: dict[RankVector, list[int]] = {}
    for index, fit in enumerate(preliminary):
        occurrences.setdefault(_rank_vector(fit.theta, threshold, caps), []).append(index)

    catalog: list[dict[str, Any]] = []
    representatives: dict[RankVector, int] = {}
    for rank, indices in occurrences.items():
        first_index = indices[0]
        interior = [
            index
            for index in indices
            if np.isfinite(max_abs(preliminary[index].theta))
            and max_abs(preliminary[index].theta) < coefficient_bound
        ]
        eligible = interior or indices
        best_index = min(
            eligible,
            key=lambda index: (preliminary[index].objective, index),
        )
        # Earlier/high-penalty singular spaces are preferred for the lowest ranks;
        # other ranks use their best-objective occurrence.
        representatives[rank] = best_index
        catalog.append(
            {
                "rank_vector": rank,
                "first_path_index": first_index,
                "best_path_index": best_index,
                "occurrence_count": len(indices),
            }
        )

    ordered_ranks = sorted(occurrences, key=lambda rank: (sum(rank), rank))
    required: list[RankVector] = []
    zero_rank = tuple(0 for _ in caps)
    if zero_rank in occurrences:
        required.append(zero_rank)
    nonzero = [rank for rank in ordered_ranks if any(rank)]
    for rank in nonzero[:2]:
        representatives[rank] = occurrences[rank][0]
        required.append(rank)
    interior_indices = [
        index
        for index, fit in enumerate(preliminary)
        if max_abs(fit.theta) < coefficient_bound and fit.converged
    ]
    if interior_indices:
        best_interior_index = min(
            interior_indices,
            key=lambda index: (preliminary[index].objective, index),
        )
        best_interior_rank = _rank_vector(
            preliminary[best_interior_index].theta, threshold, caps
        )
        representatives[best_interior_rank] = best_interior_index
        required.append(best_interior_rank)
    if ordered_ranks:
        required.append(ordered_ranks[-1])

    selected_ranks: list[RankVector] = []
    for rank in [*required, *ordered_ranks]:
        if rank not in selected_ranks:
            selected_ranks.append(rank)
        if len(selected_ranks) >= max_routes:
            break

    routes: list[CapPilotRoute] = []
    for rank in selected_ranks:
        index = representatives[rank]
        indices = occurrences[rank]
        best_index = min(
            indices, key=lambda candidate: (preliminary[candidate].objective, candidate)
        )
        fit = preliminary[index]
        routes.append(
            CapPilotRoute(
                index,
                rank,
                fit.theta,
                fit.penalty,
                fit.objective,
                indices[0],
                best_index,
                "nuclear_path",
            )
        )

    if zero_rank not in occurrences and len(routes) < max_routes:
        routes.insert(
            0,
            CapPilotRoute(
                None,
                zero_rank,
                zeros_like(preliminary[0].theta),
                None,
                None,
                None,
                None,
                "synthetic_zero",
            ),
        )
    selected = {route.rank for route in routes}
    synthetic_low_ranks = [
        tuple(1 if index == block and caps[index] > 0 else 0 for index in range(len(caps)))
        for block in range(len(caps))
        if caps[block] > 0
    ]
    minimum_route_count = min(4, max_routes)
    for rank in synthetic_low_ranks:
        if len(routes) >= minimum_route_count:
            break
        if rank in selected:
            continue
        routes.append(
            CapPilotRoute(
                None,
                rank,
                _synthetic_cap_start(preliminary[0].theta, rank, threshold),
                None,
                None,
                None,
                None,
                "synthetic_very_low_rank",
            )
        )
        selected.add(rank)
    return routes[:max_routes], sorted(
        catalog, key=lambda item: (sum(item["rank_vector"]), item["rank_vector"])
    )


def _fit_candidate(
    y: np.ndarray,
    design: Design,
    ranks: RankVector,
    starts: tuple[Coefficients, Coefficients],
    seed: int | np.random.SeedSequence,
    fit_options: dict[str, Any],
    start_objective_stability_tol: float,
    stationarity_tolerance: float,
    start_envelope_fraction: float | None = None,
    diagnostic_context: str = "candidate",
) -> tuple[FitResult, bool, list[str], dict[str, Any]]:
    prepared_starts = [adapt_initial(start, ranks) for start in starts]
    start_preparation: list[dict[str, float | bool]] = []
    if start_envelope_fraction is not None:
        prepared = []
        for start in prepared_starts:
            rescaled, diagnostics = _rescale_cap_start(
                start,
                float(fit_options.get("coefficient_bound", 10.0)),
                start_envelope_fraction,
            )
            prepared.append(rescaled)
            start_preparation.append(diagnostics)
        prepared_starts = prepared
    first = fit_fixed_rank(
        y,
        design,
        ranks,
        initial=prepared_starts[0],
        seed=seed,
        diagnostic_context=f"{diagnostic_context}:start_1",
        **fit_options,
    )
    second = fit_fixed_rank(
        y,
        design,
        ranks,
        initial=prepared_starts[1],
        seed=seed,
        diagnostic_context=f"{diagnostic_context}:start_2",
        **fit_options,
    )
    options = [first, second]
    initially_valid = [
        fit for fit in options if not fit_invalid_reasons(fit, ranks, stationarity_tolerance)
    ]
    initial_stable = bool(
        len(initially_valid) >= 2
        and abs(initially_valid[0].objective - initially_valid[1].objective)
        / max(1.0, abs(min(initially_valid[0].objective, initially_valid[1].objective)))
        <= start_objective_stability_tol
    )
    if not initial_stable:
        options.append(
            fit_fixed_rank(
                y,
                design,
                ranks,
                seed=seed,
                diagnostic_context=f"{diagnostic_context}:start_3",
                **fit_options,
            )
        )
    valid = [
        fit
        for fit in options
        if not fit_invalid_reasons(fit, ranks, stationarity_tolerance)
    ]
    ordered_valid = sorted(valid, key=lambda fit: fit.objective)
    stable = False
    best_two_gap = float("nan")
    if len(ordered_valid) >= 2:
        best_two_gap = abs(ordered_valid[0].objective - ordered_valid[1].objective) / max(
            1.0, abs(ordered_valid[0].objective)
        )
        stable = best_two_gap <= start_objective_stability_tol
    chosen = ordered_valid[0] if ordered_valid else min(options, key=lambda fit: fit.objective)
    reasons = fit_invalid_reasons(chosen, ranks, stationarity_tolerance)
    if not stable:
        reasons.append("objective_stability_failed")
    diagnostics = {
        "start_objectives": [fit.objective for fit in options],
        "start_stationarity_residuals": [fit.stationarity_residual for fit in options],
        "start_valid": [
            not fit_invalid_reasons(fit, ranks, stationarity_tolerance) for fit in options
        ],
        "best_two_objective_gap": best_two_gap,
        "objective_stability_pass": stable,
        "third_start_used": len(options) == 3,
        "start_preparation": start_preparation,
        "start_details": [
            {
                "start_number": index,
                "diagnostic_context": fit.diagnostics.get("diagnostic_context"),
                "requested_rank": fit.ranks,
                "numerical_rank": _numerical_rank_vector(fit.theta),
                "objective_initial": fit.objective_history[0],
                "objective_final": fit.objective,
                "iterations": fit.iterations,
                "converged": fit.converged,
                "stationarity_residual": fit.stationarity_residual,
                "stationarity_pass": bool(
                    fit.diagnostics.get("stationarity_pass", False)
                    if fit.diagnostics.get("constrained_fallback_used")
                    else np.isfinite(fit.stationarity_residual)
                    and fit.stationarity_residual <= stationarity_tolerance
                ),
                "initial_coefficient_envelope": fit.diagnostics.get(
                    "initial_coefficient_envelope"
                ),
                "final_coefficient_envelope": fit.diagnostics.get(
                    "final_coefficient_envelope"
                ),
                "coefficient_envelope_ratio": fit.max_envelope_ratio,
                "coefficient_bound_pass": fit.max_envelope_ratio <= 1.0 + 1e-8,
                "runtime_seconds": fit.diagnostics.get("runtime_seconds"),
                "invalid_reasons": fit_invalid_reasons(
                    fit, ranks, stationarity_tolerance
                ),
                "coefficient_envelope_history": fit.diagnostics.get(
                    "coefficient_envelope_history"
                ),
            }
            for index, fit in enumerate(options, start=1)
        ],
    }
    return chosen, len(options) == 3, reasons, diagnostics


def _rank_move_initial(
    y: np.ndarray,
    design: Design,
    fit: FitResult,
    target_ranks: RankVector,
) -> Coefficients:
    """Construct a decrease or normal-gradient increase initialization."""

    current = fit.theta
    matrices = [matrix.copy() for matrix in current.matrices()]
    changed = [
        index
        for index, (old, new) in enumerate(zip(fit.ranks, target_ranks, strict=True))
        if old != new
    ]
    if len(changed) != 1:
        return adapt_initial(current, target_ranks)
    index = changed[0]
    if target_ranks[index] < fit.ranks[index]:
        matrices[index] = truncated_matrix(matrices[index], target_ranks[index])
    else:
        residual = fitted_values(current, design) - y
        gradient = adjoint(residual, design)
        tangent = tangent_project(gradient, current, fit.ranks)
        normal = gradient.matrices()[index] - tangent.matrices()[index]
        u, singular, vt = np.linalg.svd(normal, full_matrices=False)
        if singular.size and singular[0] > 0.0:
            step = max(float(singular[0]) / (y.size), 1e-4)
            matrices[index] -= step * np.outer(u[:, 0], vt[0])
    return from_matrices(matrices, len(current.A), len(current.B))


def fit_rank_adaptive_cap_pilot(
    y: np.ndarray,
    design: Design,
    caps: RankVector,
    preliminary: list[NuclearFit],
    threshold: float,
    *,
    seed: int | np.random.SeedSequence,
    fit_options: dict[str, Any],
    stationarity_tolerance: float,
    start_objective_stability_tol: float,
    improvement_tolerance: float,
    removal_tolerance: float,
    max_steps: int,
    max_routes: int = 6,
    start_envelope_fraction: float = 0.8,
) -> tuple[FitResult, dict[str, Any]]:
    """Compute one rank-at-most-cap pilot by local rank-adaptive loss descent."""

    if not preliminary:
        raise RankPilotFailure("rank-cap pilot requires a nonempty nuclear path")
    coefficient_bound = float(fit_options.get("coefficient_bound", 10.0))
    routes, route_catalog = _cap_pilot_routes(
        preliminary,
        threshold,
        caps,
        max_routes=max_routes,
        coefficient_bound=coefficient_bound,
    )
    route_results: list[
        tuple[CandidateFit, list[dict[str, Any]], CapPilotRoute]
    ] = []
    route_attempts: list[dict[str, Any]] = []
    fit_cache: dict[
        tuple[RankVector, str, str],
        tuple[FitResult, bool, list[str], dict[str, Any]],
    ] = {}
    cache_hits = 0
    cache_misses = 0

    def cached_fit(
        ranks: RankVector,
        starts: tuple[Coefficients, Coefficients],
        diagnostic_context: str,
    ) -> tuple[FitResult, bool, list[str], dict[str, Any]]:
        nonlocal cache_hits, cache_misses
        prepared = []
        for start in starts:
            adapted = adapt_initial(start, ranks)
            rescaled, _ = _rescale_cap_start(
                adapted, coefficient_bound, start_envelope_fraction
            )
            prepared.append(rescaled)
        key = (ranks, _theta_fingerprint(prepared[0]), _theta_fingerprint(prepared[1]))
        if key in fit_cache:
            cache_hits += 1
            return fit_cache[key]
        cache_misses += 1
        result = _fit_candidate(
            y,
            design,
            ranks,
            starts,
            seed,
            fit_options,
            start_objective_stability_tol,
            stationarity_tolerance,
            start_envelope_fraction,
            diagnostic_context,
        )
        fit_cache[key] = result
        return result

    for route_number, route in enumerate(routes):
        start_theta = route.theta
        start_rank = route.rank
        if route.path_index is None:
            second_theta = zeros_like(start_theta)
            second_path_index = None
        else:
            alternatives = [
                index
                for index in (route.first_path_index, route.best_path_index)
                if index is not None and index != route.path_index
            ]
            second_path_index = (
                alternatives[0]
                if alternatives
                else max(0, route.path_index - 1)
            )
            second_theta = preliminary[second_path_index].theta
        prepared_start, preparation = _rescale_cap_start(
            adapt_initial(start_theta, start_rank),
            coefficient_bound,
            start_envelope_fraction,
        )
        fit, third, reasons, start_diagnostics = cached_fit(
            start_rank,
            (start_theta, second_theta),
            f"cap_pilot_route_{route_number}:initial_rank_{start_rank}",
        )
        current = CandidateFit(
            start_rank,
            fit,
            float("nan"),
            model_dimension(start_rank, *y.shape),
            [
                f"nuclear_path_{route.path_index}"
                if route.path_index is not None
                else "synthetic_zero_start"
            ],
            third,
            not reasons,
            reasons,
            start_diagnostics,
        )
        path = [
            {
                "ranks": start_rank,
                "objective": fit.objective,
                "valid": current.valid,
                "move": "start",
            }
        ]
        collapsed_rank = _numerical_rank_vector(fit.theta)
        collapse_only = bool(reasons) and all(
            reason.startswith("numerical_rank_support:")
            or reason == "objective_stability_failed"
            for reason in reasons
        )
        if collapse_only and collapsed_rank != start_rank:
            collapsed_fit, collapsed_third, collapsed_reasons, collapsed_starts = cached_fit(
                collapsed_rank,
                (fit.theta, start_theta),
                f"cap_pilot_route_{route_number}:collapse_rank_{collapsed_rank}",
            )
            current = CandidateFit(
                collapsed_rank,
                collapsed_fit,
                float("nan"),
                model_dimension(collapsed_rank, *y.shape),
                [f"numerical_collapse_of_{start_rank}"],
                collapsed_third,
                not collapsed_reasons,
                collapsed_reasons,
                collapsed_starts,
            )
            path.append(
                {
                    "ranks": collapsed_rank,
                    "objective": collapsed_fit.objective,
                    "valid": current.valid,
                    "move": "numerical_rank_collapse",
                }
            )
        route_attempt = {
            "route_type": "original",
            "route_number": route_number,
            "route_source": route.source,
            "nuclear_path_index": route.path_index,
            "secondary_path_index": second_path_index,
            "nuclear_penalty": route.penalty,
            "nuclear_objective": route.nuclear_objective,
            "first_path_index_for_rank": route.first_path_index,
            "best_path_index_for_rank": route.best_path_index,
            "start_rank_vector": start_rank,
            "original_start_max_abs_coefficient": preparation[
                "original_max_abs_coefficient"
            ],
            "start_common_scale": preparation["applied_common_scale"],
            "rescaled_start_max_abs_coefficient": preparation[
                "rescaled_max_abs_coefficient"
            ],
            "start_common_rescaling_used": preparation["common_rescaling_used"],
            "prepared_start_fingerprint": _theta_fingerprint(prepared_start),
            "start_fit_reasons": reasons,
            "start_fit_diagnostics": start_diagnostics,
            "path": path,
        }
        route_attempts.append(route_attempt)
        visited = {start_rank, current.ranks}
        for _ in range(max_steps):
            if not current.valid:
                break
            neighbors: list[CandidateFit] = []
            for ranks in sorted(one_coordinate_neighbors(current.ranks, caps) - visited):
                move_initial = _rank_move_initial(y, design, current.fit, ranks)
                closest = _closest_preliminary(ranks, preliminary, caps, threshold)
                neighbor_fit, used_third, neighbor_reasons, neighbor_starts = cached_fit(
                    ranks,
                    (move_initial, closest),
                    f"cap_pilot_route_{route_number}:neighbor_rank_{ranks}",
                )
                neighbors.append(
                    CandidateFit(
                        ranks,
                        neighbor_fit,
                        float("nan"),
                        model_dimension(ranks, *y.shape),
                        [f"adaptive_neighbor_of_{current.ranks}"],
                        used_third,
                        not neighbor_reasons,
                        neighbor_reasons,
                        neighbor_starts,
                    )
                )
                visited.add(ranks)
            valid_neighbors = [item for item in neighbors if item.valid]
            objective_scale = max(1.0, abs(current.fit.objective))
            improvements = [
                item
                for item in valid_neighbors
                if item.fit.objective
                < current.fit.objective - improvement_tolerance * objective_scale
            ]
            move = "loss_improvement"
            if improvements:
                chosen = min(improvements, key=lambda item: (item.fit.objective, item.dimension))
            else:
                removals = [
                    item
                    for item in valid_neighbors
                    if item.dimension < current.dimension
                    and item.fit.objective
                    <= current.fit.objective + removal_tolerance * objective_scale
                ]
                if not removals:
                    break
                chosen = min(removals, key=lambda item: (item.dimension, item.fit.objective))
                move = "numerically_redundant_removal"
            current = chosen
            path.append(
                {
                    "ranks": current.ranks,
                    "objective": current.fit.objective,
                    "valid": current.valid,
                    "move": move,
                }
            )
        if current.valid:
            route_results.append((current, path, route))
        route_attempt.update(
            {
                "final_rank_vector": current.ranks,
                "final_numerical_rank_vector": _numerical_rank_vector(current.fit.theta),
                "final_thresholded_rank_vector": _rank_vector(
                    current.fit.theta, threshold, caps
                ),
                "final_objective": current.fit.objective,
                "final_stationarity_residual": current.fit.stationarity_residual,
                "final_max_envelope_ratio": current.fit.max_envelope_ratio,
                "final_converged": current.fit.converged,
                "final_valid": current.valid,
                "final_reasons": current.invalid_reasons,
            }
        )

    ordered = sorted(route_results, key=lambda item: item[0].fit.objective)
    best_objective = ordered[0][0].fit.objective if ordered else float("nan")
    second_objective = ordered[1][0].fit.objective if len(ordered) >= 2 else float("nan")
    outer_gap = (
        abs(best_objective - second_objective) / max(1.0, abs(best_objective))
        if len(ordered) >= 2
        else float("nan")
    )
    original_stable = len(ordered) >= 2 and outer_gap <= start_objective_stability_tol
    confirmation_fits: list[FitResult] = []
    confirmation_records: list[dict[str, Any]] = []
    confirmation_pass = False
    acceptance_basis = (
        "best_valid_route_agreement"
        if original_stable
        else "best_valid_route_disagreement_warning"
    )
    if ordered and not original_stable:
        confirmation_fits, confirmation_records, confirmation_pass = _confirm_best_basin(
            y,
            design,
            ordered[0][0].fit,
            ordered[0][0].ranks,
            fit_options=fit_options,
            seed=seed,
            stationarity_tolerance=stationarity_tolerance,
            objective_tolerance=start_objective_stability_tol,
            start_envelope_fraction=start_envelope_fraction,
            diagnostic_context="cap_pilot_best_basin",
        )
    confirmation_valid = [record for record in confirmation_records if record["valid"]]
    confirmation_matching = [
        record for record in confirmation_records if record["confirmation_pass"]
    ]
    multistart_agreement = original_stable or confirmation_pass
    if confirmation_pass and not original_stable:
        acceptance_basis = "best_valid_route_confirmed_agreement"
    common_diagnostics = {
        "attempted_route_count": len(routes),
        "valid_route_count": len(ordered),
        "stable_route_count": (
            sum(
                abs(item[0].fit.objective - best_objective)
                / max(1.0, abs(best_objective))
                <= start_objective_stability_tol
                for item in ordered
            )
            if ordered
            else 0
        ),
        "best_two_objective_gap": outer_gap,
        "objective_stability_pass": multistart_agreement,
        "best_valid_objective": best_objective,
        "second_best_valid_objective": second_objective,
        "normalized_objective_gap": outer_gap,
        "multistart_objective_agreement": multistart_agreement,
        "pilot_multistart_disagreement": not multistart_agreement,
        "basin_confirmation_attempted": bool(ordered and not original_stable),
        "basin_confirmation_success": confirmation_pass,
        "outer_start_attempts": route_attempts,
        "basin_confirmation_attempts": confirmation_records,
        "original_best_objective": best_objective,
        "original_second_best_objective": second_objective,
        "original_stability_gap": outer_gap,
        "confirmation_best_objective": min(
            (record["final_objective"] for record in confirmation_valid),
            default=float("nan"),
        ),
        "number_confirmation_valid": len(confirmation_valid),
        "number_confirmation_matching_best": len(confirmation_matching),
        "final_pilot_acceptance_basis": acceptance_basis,
    }
    if not ordered:
        common_diagnostics["final_pilot_acceptance_basis"] = "failure"
        raise RankPilotFailure(
            "rank-adaptive cap pilot has no valid outer route", common_diagnostics
        )
    chosen = ordered[0][0]
    stable_routes = [
        item
        for item in ordered
        if abs(item[0].fit.objective - best_objective)
        / max(1.0, abs(best_objective))
        <= start_objective_stability_tol
    ]
    stable_numerical_ranks = [
        _numerical_rank_vector(item[0].fit.theta) for item in stable_routes
    ]
    stable_thresholded_ranks = [
        _rank_vector(item[0].fit.theta, threshold, caps) for item in stable_routes
    ]
    return chosen.fit, {
        "algorithm": "rank_adaptive_at_most_cap",
        "route_catalog": route_catalog,
        "attempted_route_count": len(routes),
        "valid_route_count": len(ordered),
        "stable_route_count": len(stable_routes),
        "outer_start_rank_vectors": [item[2].rank for item in ordered],
        "outer_start_path_indices": [item[2].path_index for item in ordered],
        "outer_final_rank_vectors": [item[0].ranks for item in ordered],
        "outer_objectives": [item[0].fit.objective for item in ordered],
        "outer_paths": [item[1] for item in ordered],
        "outer_start_attempts": route_attempts,
        "basin_confirmation_attempts": confirmation_records,
        "best_two_objective_gap": outer_gap,
        "objective_stability_pass": multistart_agreement,
        "best_valid_objective": best_objective,
        "second_best_valid_objective": second_objective,
        "normalized_objective_gap": outer_gap,
        "multistart_objective_agreement": multistart_agreement,
        "pilot_multistart_disagreement": not multistart_agreement,
        "basin_confirmation_attempted": not original_stable,
        "basin_confirmation_success": confirmation_pass,
        "original_best_objective": best_objective,
        "original_second_best_objective": second_objective,
        "original_stability_gap": outer_gap,
        "confirmation_best_objective": common_diagnostics[
            "confirmation_best_objective"
        ],
        "number_confirmation_valid": len(confirmation_valid),
        "number_confirmation_matching_best": len(confirmation_matching),
        "final_pilot_acceptance_basis": acceptance_basis,
        "stable_final_numerical_rank_vectors": stable_numerical_ranks,
        "stable_final_thresholded_rank_vectors": stable_thresholded_ranks,
        "stable_final_numerical_ranks_agree": len(set(stable_numerical_ranks)) == 1,
        "stable_final_thresholded_ranks_agree": len(set(stable_thresholded_ranks)) == 1,
        "route_fit_cache_hits": cache_hits,
        "route_fit_cache_misses": cache_misses,
        "numerical_rank_before_thresholding": _numerical_rank_vector(chosen.fit.theta),
        "objective": chosen.fit.objective,
        "stationarity_residual": chosen.fit.stationarity_residual,
        "max_envelope_ratio": chosen.fit.max_envelope_ratio,
        "starting_values": (
            "up to six distinct-rank deterministic nuclear-path routes, including "
            "low-rank and zero routes when admissible"
        ),
        "start_envelope_fraction": start_envelope_fraction,
        "coefficient_bound": coefficient_bound,
    }


def fit_revision10_spectral_pilot(
    y: np.ndarray,
    design: Design,
    reporting_caps: RankVector,
    *,
    seed: int | np.random.SeedSequence,
    fit_options: dict[str, Any],
    stationarity_tolerance: float,
    start_objective_stability_tol: float,
) -> tuple[FitResult, dict[str, Any]]:
    """Fit one joint at-most-cap+1 pilot with maintained numerical diagnostics."""

    pilot_caps = tuple(cap + 1 for cap in reporting_caps)
    if any(cap < 0 for cap in reporting_caps) or any(
        cap > min(y.shape) for cap in pilot_caps
    ):
        raise ValueError("reporting caps must admit their cap+1 pilot ranks")
    rng = np.random.default_rng(seed)
    additional = [int(value) for value in rng.integers(0, 2**32 - 1, size=2)]
    start_seeds: list[int | np.random.SeedSequence] = [seed, *additional]
    fits: list[FitResult] = []
    records: list[dict[str, Any]] = []
    for index, start_seed in enumerate(start_seeds, start=1):
        fit = fit_fixed_rank(
            y,
            design,
            pilot_caps,
            seed=start_seed,
            diagnostic_context=f"revision10_spectral_pilot:start_{index}",
            **fit_options,
        )
        reasons = fit_invalid_reasons(
            fit,
            pilot_caps,
            stationarity_tolerance,
            require_exact_numerical_rank=False,
        )
        fits.append(fit)
        records.append(
            {
                "start_id": f"deterministic_start_{index}",
                "seed": (
                    f"seed_sequence:{start_seed.entropy}:{start_seed.spawn_key}"
                    if isinstance(start_seed, np.random.SeedSequence)
                    else str(start_seed)
                ),
                "objective": fit.objective,
                "finite_objective": bool(np.isfinite(fit.objective)),
                "feasible": "constrained_feasibility_failure" not in reasons,
                "converged": fit.converged,
                "stationarity_residual": fit.stationarity_residual,
                "stationarity_pass": "stationarity_high" not in reasons
                and "constrained_optimality_failure" not in reasons,
                "max_envelope_ratio": fit.max_envelope_ratio,
                "boundary_active": bool(fit.diagnostics.get("boundary_active", False)),
                "termination_status": fit.diagnostics.get(
                    "constrained_solver_status", "unknown"
                ),
                "numerical_rank_vector": _numerical_rank_vector(fit.theta),
                "valid": not reasons,
                "invalid_reasons": reasons,
            }
        )
    valid = sorted(
        (fit for fit, record in zip(fits, records, strict=True) if record["valid"]),
        key=lambda fit: fit.objective,
    )
    best = valid[0] if valid else min(fits, key=lambda fit: fit.objective)
    second_objective = valid[1].objective if len(valid) >= 2 else float("nan")
    objective_gap = (
        abs(second_objective - best.objective) / max(1.0, abs(best.objective))
        if len(valid) >= 2
        else float("nan")
    )
    stable = len(valid) >= 2 and objective_gap <= start_objective_stability_tol
    diagnostics = {
        "algorithm": "joint_at_most_cap_plus_one_three_deterministic_starts",
        "reporting_rank_caps": reporting_caps,
        "pilot_rank_caps": pilot_caps,
        "start_objective_stability_tolerance": start_objective_stability_tol,
        "maintained_start_count": 3,
        "third_start_used": True,
        "start_records": records,
        "all_start_objectives": [record["objective"] for record in records],
        "all_start_stationarity_residuals": [
            record["stationarity_residual"] for record in records
        ],
        "valid_start_count": len(valid),
        "best_objective": best.objective,
        "second_best_valid_objective": second_objective,
        "best_two_objective_gap": objective_gap,
        "objective_stability_pass": stable,
        "feasibility": bool(
            next(
                record["feasible"]
                for fit, record in zip(fits, records, strict=True)
                if fit is best
            )
        ),
        "finite_objective": bool(np.isfinite(best.objective)),
        "stationarity_residual": best.stationarity_residual,
        "coefficient_box_activity": bool(
            best.diagnostics.get("boundary_active", False)
        ),
        "max_envelope_ratio": best.max_envelope_ratio,
        "termination_status": best.diagnostics.get("constrained_solver_status", "unknown"),
        "numerical_rank_vector": _numerical_rank_vector(best.theta),
        "starting_values": [record["start_id"] for record in records],
        "global_objective_gap_certified": False,
    }
    if not stable:
        raise RankPilotFailure(
            "rank_selection_numerically_unresolved: cap+1 pilot failed maintained acceptance diagnostics",
            diagnostics,
        )
    return best, diagnostics


def select_ranks(
    y: np.ndarray,
    design: Design,
    caps: RankVector,
    *,
    seed: int | np.random.SeedSequence = 0,
    coefficient_bound: float = 10.0,
    max_sweeps: int = 200,
    objective_rtol: float = 1e-8,
    stationarity_tol: float = 1e-6,
    lstsq_rcond: float = 1e-10,
    interior_numerical_tolerance: float = 1e-8,
    constraint_tolerance: float = 1e-8,
    constrained_kkt_tolerance: float = 1e-4,
    constrained_subproblem_tolerance: float = 1e-10,
    constrained_subproblem_max_iterations: int = 200,
    start_objective_stability_tol: float = 1e-6,
    cap_pilot_start_envelope_fraction: float = 0.8,
) -> Revision10RankSelectionResult:
    """Apply the frozen blockwise ridge-ratio selector in Revision-10 Section 4.5."""

    n, t = y.shape
    if design.shape != y.shape:
        raise ValueError("response and design shapes differ")
    if len(caps) != len(design.y_lags) + len(design.x) + 1:
        raise ValueError("rank caps do not match coefficient blocks")
    fit_options = {
        "coefficient_bound": coefficient_bound,
        "max_sweeps": max_sweeps,
        "objective_rtol": objective_rtol,
        "stationarity_tol": stationarity_tol,
        "lstsq_rcond": lstsq_rcond,
        "interior_numerical_tolerance": interior_numerical_tolerance,
        "constraint_tolerance": constraint_tolerance,
        "constrained_kkt_tolerance": constrained_kkt_tolerance,
        "constrained_subproblem_tolerance": constrained_subproblem_tolerance,
        "constrained_subproblem_max_iterations": constrained_subproblem_max_iterations,
    }
    pilot, pilot_diagnostics = fit_revision10_spectral_pilot(
        y,
        design,
        caps,
        seed=seed,
        fit_options=fit_options,
        stationarity_tolerance=stationarity_tol,
        start_objective_stability_tol=start_objective_stability_tol,
    )
    names = revision10_block_names(design)
    weights = revision10_scale_weights(design)
    reference = weights[0]
    block_records: dict[str, dict[str, Any]] = {}
    all_ratios: list[tuple[float, ...]] = []
    for name, matrix, weight, cap in zip(
        names, pilot.theta.matrices(), weights, caps, strict=True
    ):
        singular, normalized = revision10_normalized_spectrum(
            matrix,
            block_weight=weight,
            reference_weight=reference,
            n=n,
            t=t,
            count=cap + 1,
        )
        ratios = revision10_ridge_ratios(
            normalized, reporting_cap=cap, n=n, t=t
        )
        selected_rank = revision10_select_block_rank(ratios)
        ordered_ratios = sorted(ratios)
        second = ordered_ratios[1] if len(ordered_ratios) >= 2 else None
        block_records[name] = {
            "weight": weight,
            "pilot_singular_values_through_cap_plus_one": singular,
            "normalized_lambda_hat_through_cap_plus_one": normalized,
            "ratios_R_M_0_through_cap": ratios,
            "minimum_ratio": ordered_ratios[0],
            "second_smallest_ratio": second,
            "ratio_gap": second - ordered_ratios[0] if second is not None else None,
            "selected_rank": selected_rank,
        }
        all_ratios.append(ratios)
    selected_ranks = revision10_assemble_rank_vector(all_ratios)
    final_fit, final_diagnostics = fit_fixed_rank_multistart(
        y,
        design,
        selected_ranks,
        seed=seed,
        fit_options=fit_options,
        stationarity_tolerance=stationarity_tol,
        start_objective_stability_tol=start_objective_stability_tol,
        start_envelope_fraction=cap_pilot_start_envelope_fraction,
    )
    final_reasons = fit_invalid_reasons(final_fit, selected_ranks, stationarity_tol)
    if not final_diagnostics["objective_stability_pass"] or final_reasons:
        raise FinalPostRefitFailure(
            "selected_rank_post_refit_numerically_unresolved",
            {
                "selected_rank_vector": selected_ranks,
                "invalid_reasons": final_reasons,
                "multistart": final_diagnostics,
            },
        )
    pilot_diagnostics["singular_values_through_cap_plus_one"] = {
        name: block_records[name]["pilot_singular_values_through_cap_plus_one"]
        for name in names
    }
    diagnostics = {
        "rank_selector_method": "revision10_ridge_ratio",
        "reporting_rank_caps": caps,
        "pilot_rank_caps": tuple(cap + 1 for cap in caps),
        "block_order": names,
        "scale_weights": dict(zip(names, weights, strict=True)),
        "reference_weight_w_A1": reference,
        "a_NT": revision10_ridge(n, t),
        "blocks": block_records,
        "selected_rank_by_block": {
            name: block_records[name]["selected_rank"] for name in names
        },
        "selected_rank_vector": selected_ranks,
        "selected_rank_at_reporting_cap": any(
            rank == cap for rank, cap in zip(selected_ranks, caps, strict=True)
        ),
        "pilot": pilot_diagnostics,
        "rank_selection_numerically_unresolved": False,
        "final_selected_rank_post_refit_status": "success",
        "final_post_refit_multistart": final_diagnostics,
    }
    return Revision10RankSelectionResult(selected_ranks, pilot, final_fit, diagnostics)


def select_ranks_revision9(
    y: np.ndarray,
    design: Design,
    caps: RankVector,
    *,
    seed: int | np.random.SeedSequence = 0,
    coefficient_bound: float = 10.0,
    max_sweeps: int = 200,
    objective_rtol: float = 1e-8,
    stationarity_tol: float = 1e-6,
    lstsq_rcond: float = 1e-10,
    interior_numerical_tolerance: float = 1e-8,
    constraint_tolerance: float = 1e-8,
    constrained_kkt_tolerance: float = 1e-4,
    constrained_subproblem_tolerance: float = 1e-10,
    constrained_subproblem_max_iterations: int = 200,
    nuclear_gamma: float = 0.8,
    nuclear_epsilon: float = 0.01,
    nuclear_max_iter: int = 500,
    nuclear_tol: float = 1e-7,
    dykstra_max_iter: int = 100,
    dykstra_tol: float = 1e-9,
    eta_for_penalty: float = 4.0,
    spatial_dimension: int = 1,
    ic_multiplier: float = 1.0,
    threshold_multiplier: float = 1.0,
    start_objective_stability_tol: float = 1e-6,
    rank_adaptive_improvement_tol: float = 1e-7,
    rank_adaptive_removal_tol: float = 1e-7,
    rank_adaptive_max_steps: int = 12,
    rank_adaptive_max_routes: int = 6,
    cap_pilot_start_envelope_fraction: float = 0.8,
    dense_nuclear_gamma: float = float(np.sqrt(0.8)),
    threshold_sensitivity_multipliers: list[float] | tuple[float, ...] = (0.5, 2.0),
    ic_sensitivity_multipliers: list[float] | tuple[float, ...] = (0.5, 2.0),
    larger_rank_caps: list[int] | tuple[int, ...] = (4, 4, 4),
    compute_rank_sensitivities: bool = True,
    compute_dense_grid_sensitivity: bool = True,
    compute_larger_cap_sensitivity: bool = True,
) -> RankSelectionResult:
    """Preserved historical Revision-9 IC selector; never the primary Revision-10 path."""

    n, t = y.shape
    if len(caps) != len(design.y_lags) + len(design.x) + 1:
        raise ValueError("rank caps do not match coefficient blocks")
    nuclear_options = {
        "coefficient_bound": coefficient_bound,
        "max_iter": nuclear_max_iter,
        "tolerance": nuclear_tol,
        "dykstra_max_iter": dykstra_max_iter,
        "dykstra_tolerance": dykstra_tol,
    }
    preliminary = nuclear_path(
        y, design, gamma=nuclear_gamma, epsilon=nuclear_epsilon, **nuclear_options
    )
    fit_options = {
        "coefficient_bound": coefficient_bound,
        "max_sweeps": max_sweeps,
        "objective_rtol": objective_rtol,
        "stationarity_tol": stationarity_tol,
        "lstsq_rcond": lstsq_rcond,
        "interior_numerical_tolerance": interior_numerical_tolerance,
        "constraint_tolerance": constraint_tolerance,
        "constrained_kkt_tolerance": constrained_kkt_tolerance,
        "constrained_subproblem_tolerance": constrained_subproblem_tolerance,
        "constrained_subproblem_max_iterations": constrained_subproblem_max_iterations,
    }
    baseline_threshold = threshold_multiplier * np.sqrt(n * t) / np.log(n * t)
    cap_fit, cap_pilot_diagnostics = fit_rank_adaptive_cap_pilot(
        y,
        design,
        caps,
        preliminary,
        baseline_threshold,
        seed=seed,
        fit_options=fit_options,
        stationarity_tolerance=stationarity_tol,
        start_objective_stability_tol=start_objective_stability_tol,
        improvement_tolerance=rank_adaptive_improvement_tol,
        removal_tolerance=rank_adaptive_removal_tol,
        max_steps=rank_adaptive_max_steps,
        max_routes=rank_adaptive_max_routes,
        start_envelope_fraction=cap_pilot_start_envelope_fraction,
    )
    candidate_ranks, sources, path_ranks = build_candidates(
        preliminary, cap_fit, caps, baseline_threshold
    )
    initial_count = len(candidate_ranks)
    evaluated: dict[RankVector, CandidateFit] = {}

    def evaluate(
        ranks: RankVector,
        *,
        source_labels: list[str] | None = None,
        alternative_cap: Coefficients | None = None,
        preliminary_path: list[NuclearFit] | None = None,
        candidate_caps: RankVector | None = None,
        candidate_threshold: float | None = None,
    ) -> CandidateFit:
        if ranks in evaluated:
            if source_labels:
                evaluated[ranks].sources = sorted(set(evaluated[ranks].sources + source_labels))
            return evaluated[ranks]
        path = preliminary_path or preliminary
        local_caps = candidate_caps or caps
        local_threshold = (
            baseline_threshold if candidate_threshold is None else candidate_threshold
        )
        closest = _closest_preliminary(ranks, path, local_caps, local_threshold)
        fit, third, reasons, start_diagnostics = _fit_candidate(
            y,
            design,
            ranks,
            (closest, alternative_cap or cap_fit.theta),
            seed,
            fit_options,
            start_objective_stability_tol,
            stationarity_tol,
            diagnostic_context=f"post_refit_rank_{ranks}",
        )
        dimension = model_dimension(ranks, n, t)
        valid = not reasons
        ic = (
            information_criterion(
                fit.objective,
                ranks,
                n,
                t,
                eta_for_penalty,
                spatial_dimension,
                ic_multiplier,
            )
            if valid
            else float("inf")
        )
        result = CandidateFit(
            ranks,
            fit,
            ic,
            dimension,
            source_labels or sources.get(ranks, ["local_completion"]),
            third,
            valid,
            reasons,
            start_diagnostics,
        )
        evaluated[ranks] = result
        return result

    for ranks in sorted(candidate_ranks):
        evaluate(ranks)

    def best(multiplier: float = ic_multiplier, pool: set[RankVector] | None = None) -> CandidateFit:
        eligible = [
            item
            for ranks, item in evaluated.items()
            if item.valid and (pool is None or ranks in pool)
        ]
        if not eligible:
            raise RankSelectionFailure("no valid candidate post-refit remains")
        return min(
            eligible,
            key=lambda item: (
                information_criterion(
                    item.fit.objective,
                    item.ranks,
                    n,
                    t,
                    eta_for_penalty,
                    spatial_dimension,
                    multiplier,
                ),
                item.dimension,
                item.ranks,
            ),
        )

    def complete_selection(
        initial_pool: set[RankVector],
        local_caps: RankVector,
        *,
        multiplier: float,
        label: str,
        alternative_cap: Coefficients | None = None,
        preliminary_path: list[NuclearFit] | None = None,
        threshold: float = baseline_threshold,
    ) -> tuple[CandidateFit, set[RankVector]]:
        pool = set(initial_pool)
        for ranks in sorted(pool):
            evaluate(
                ranks,
                source_labels=[label],
                alternative_cap=alternative_cap,
                preliminary_path=preliminary_path,
                candidate_caps=local_caps,
                candidate_threshold=threshold,
            )
        while True:
            selected_local = best(multiplier=multiplier, pool=pool)
            missing = one_coordinate_neighbors(selected_local.ranks, local_caps) - pool
            for ranks in sorted(missing):
                evaluate(
                    ranks,
                    source_labels=[f"{label}_neighbor_of_{selected_local.ranks}"],
                    alternative_cap=alternative_cap,
                    preliminary_path=preliminary_path,
                    candidate_caps=local_caps,
                    candidate_threshold=threshold,
                )
            pool.update(missing)
            improved = best(multiplier=multiplier, pool=pool)
            if improved.ranks == selected_local.ranks:
                return improved, pool

    selected, baseline_candidate_vectors = complete_selection(
        set(candidate_ranks),
        caps,
        multiplier=ic_multiplier,
        label="baseline_local_completion",
        alternative_cap=cap_fit.theta,
    )
    baseline_final_count = len(evaluated)
    baseline_sources = {
        ranks: list(item.sources) for ranks, item in evaluated.items()
    }
    selected = best()
    ordered = sorted(
        (item for item in evaluated.values() if item.valid),
        key=lambda item: (item.ic, item.dimension, item.ranks),
    )
    selected_neighbors = [
        evaluated[rank]
        for rank in one_coordinate_neighbors(selected.ranks, caps)
        if rank in evaluated and evaluated[rank].valid
    ]
    neighbor_gaps = {
        str(item.ranks): float(item.ic - selected.ic) for item in selected_neighbors
    }
    best_neighbor = min(selected_neighbors, key=lambda item: item.ic) if selected_neighbors else None

    sensitivity: dict[str, Any] = {}
    if compute_rank_sensitivities:
        threshold_results = {}
        for multiplier in threshold_sensitivity_multipliers:
            threshold = float(multiplier) * np.sqrt(n * t) / np.log(n * t)
            threshold_cap_fit, threshold_cap_diagnostics = fit_rank_adaptive_cap_pilot(
                y,
                design,
                caps,
                preliminary,
                threshold,
                seed=seed,
                fit_options=fit_options,
                stationarity_tolerance=stationarity_tol,
                start_objective_stability_tol=start_objective_stability_tol,
                improvement_tolerance=rank_adaptive_improvement_tol,
                removal_tolerance=rank_adaptive_removal_tol,
                max_steps=rank_adaptive_max_steps,
                max_routes=rank_adaptive_max_routes,
                start_envelope_fraction=cap_pilot_start_envelope_fraction,
            )
            ranks_set, _, proposals = build_candidates(
                preliminary, threshold_cap_fit, caps, threshold
            )
            chosen, completed_pool = complete_selection(
                ranks_set,
                caps,
                multiplier=ic_multiplier,
                label=f"threshold_{multiplier}",
                alternative_cap=threshold_cap_fit.theta,
                threshold=threshold,
            )
            threshold_results[str(multiplier)] = {
                "selected_rank": chosen.ranks,
                "candidate_count_initial": len(ranks_set),
                "candidate_count_final": len(completed_pool),
                "candidate_rank_vectors": sorted(completed_pool),
                "path_proposals": proposals,
                "local_completion_applied": True,
                "cap_pilot_algorithm": threshold_cap_diagnostics["algorithm"],
            }
        penalty_results = {}
        for multiplier in ic_sensitivity_multipliers:
            chosen, completed_pool = complete_selection(
                set(baseline_candidate_vectors),
                caps,
                multiplier=ic_multiplier * float(multiplier),
                label=f"ic_{multiplier}",
                alternative_cap=cap_fit.theta,
            )
            penalty_results[str(multiplier)] = {
                "selected_rank": chosen.ranks,
                "candidate_count_final": len(completed_pool),
                "candidate_rank_vectors": sorted(completed_pool),
                "local_completion_applied": True,
            }
        dense_ranks: set[RankVector] = set()
        dense_proposals: list[RankVector] = []
        dense_converged: list[bool] = []
        dense_selected: RankVector | None = None
        dense_completed_count = 0
        dense_cap_algorithm: str | None = None
        if compute_dense_grid_sensitivity:
            dense = nuclear_path(
                y,
                design,
                gamma=dense_nuclear_gamma,
                epsilon=nuclear_epsilon,
                **nuclear_options,
            )
            dense_cap_fit, dense_cap_diagnostics = fit_rank_adaptive_cap_pilot(
                y,
                design,
                caps,
                dense,
                baseline_threshold,
                seed=seed,
                fit_options=fit_options,
                stationarity_tolerance=stationarity_tol,
                start_objective_stability_tol=start_objective_stability_tol,
                improvement_tolerance=rank_adaptive_improvement_tol,
                removal_tolerance=rank_adaptive_removal_tol,
                max_steps=rank_adaptive_max_steps,
                max_routes=rank_adaptive_max_routes,
                start_envelope_fraction=cap_pilot_start_envelope_fraction,
            )
            dense_ranks, _, dense_proposals = build_candidates(
                dense, dense_cap_fit, caps, baseline_threshold
            )
            dense_choice, dense_completed = complete_selection(
                dense_ranks,
                caps,
                multiplier=ic_multiplier,
                label="dense_grid",
                alternative_cap=dense_cap_fit.theta,
                preliminary_path=dense,
            )
            dense_selected = dense_choice.ranks
            dense_completed_count = len(dense_completed)
            dense_cap_algorithm = dense_cap_diagnostics["algorithm"]
            dense_converged = [fit.converged for fit in dense]
        larger_caps_tuple = tuple(int(value) for value in larger_rank_caps)
        larger_valid = False
        larger_selected: RankVector | None = None
        larger_completed_count = 0
        larger_cap_algorithm: str | None = None
        if compute_larger_cap_sensitivity:
            larger_cap_fit, larger_cap_diagnostics = fit_rank_adaptive_cap_pilot(
                y,
                design,
                larger_caps_tuple,
                preliminary,
                baseline_threshold,
                seed=seed,
                fit_options=fit_options,
                stationarity_tolerance=stationarity_tol,
                start_objective_stability_tol=start_objective_stability_tol,
                improvement_tolerance=rank_adaptive_improvement_tol,
                removal_tolerance=rank_adaptive_removal_tol,
                max_steps=rank_adaptive_max_steps,
                max_routes=rank_adaptive_max_routes,
                start_envelope_fraction=cap_pilot_start_envelope_fraction,
            )
            larger_valid = True
            larger_cap_algorithm = larger_cap_diagnostics["algorithm"]
            if larger_valid:
                larger_set, _, _ = build_candidates(
                    preliminary, larger_cap_fit, larger_caps_tuple, baseline_threshold
                )
                larger_choice, larger_completed = complete_selection(
                    larger_set,
                    larger_caps_tuple,
                    multiplier=ic_multiplier,
                    label="larger_cap",
                    alternative_cap=larger_cap_fit.theta,
                )
                larger_selected = larger_choice.ranks
                larger_completed_count = len(larger_completed)
        sensitivity = {
            "threshold_multipliers": threshold_results,
            "ic_multipliers": penalty_results,
            "dense_grid_selected_rank": dense_selected,
            "dense_grid_path_proposals": dense_proposals,
            "dense_grid_converged": dense_converged,
            "dense_grid_candidate_count_final": dense_completed_count,
            "dense_grid_local_completion_applied": bool(dense_ranks),
            "dense_grid_cap_pilot_algorithm": dense_cap_algorithm,
            "larger_caps": larger_caps_tuple,
            "larger_cap_valid": larger_valid,
            "larger_cap_selected_rank": larger_selected,
            "larger_cap_candidate_count_final": larger_completed_count,
            "larger_cap_local_completion_applied": compute_larger_cap_sensitivity,
            "larger_cap_pilot_algorithm": larger_cap_algorithm,
        }

    all_records = sorted(evaluated.values(), key=lambda item: (not item.valid, item.ic, item.ranks))
    cap_thresholded = _rank_vector(cap_fit.theta, baseline_threshold, caps)
    diagnostics = {
        "revision8_kappa": revision8_kappa(
            n,
            t,
            eta_for_penalty=eta_for_penalty,
            spatial_dimension=spatial_dimension,
            multiplier=ic_multiplier,
        ),
        "spatial_dimension": spatial_dimension,
        "threshold": float(baseline_threshold),
        "nuclear_path_rank_proposals": path_ranks,
        "distinct_nuclear_path_rank_proposals": sorted(set(path_ranks)),
        "rank_cap_thresholded_vector": cap_thresholded,
        "candidate_rank_vectors": sorted(baseline_candidate_vectors),
        "candidate_count_initial": initial_count,
        "candidate_count_final": baseline_final_count,
        "smallest_ic": selected.ic,
        "second_smallest_ic": ordered[1].ic if len(ordered) > 1 else None,
        "ic_gap": ordered[1].ic - selected.ic if len(ordered) > 1 else None,
        "ic_gaps_to_neighbors": neighbor_gaps,
        "best_neighbor_rank": best_neighbor.ranks if best_neighbor else None,
        "selected_to_cap_margins": tuple(
            cap - rank for rank, cap in zip(selected.ranks, caps, strict=True)
        ),
        "selected_rank_at_cap": any(
            rank == cap for rank, cap in zip(selected.ranks, caps, strict=True)
        ),
        "third_start_events": sum(item.third_start_used for item in evaluated.values()),
        "nuclear_path_converged": [fit.converged for fit in preliminary],
        "nuclear_path_iterations": [fit.iterations for fit in preliminary],
        "candidate_records": [
            {
                "ranks": item.ranks,
                "objective": item.fit.objective,
                "ic": item.ic,
                "dimension": item.dimension,
                "converged": item.fit.converged,
                "stationarity_residual": item.fit.stationarity_residual,
                "max_envelope_ratio": item.fit.max_envelope_ratio,
                "numerical_rank_vector": _numerical_rank_vector(item.fit.theta),
                "valid": item.valid,
                "invalid_reasons": item.invalid_reasons,
                "sources": item.sources,
                "third_start_used": item.third_start_used,
                **item.start_diagnostics,
            }
            for item in all_records
        ],
        "baseline_candidate_records": [
            {
                "ranks": item.ranks,
                "ic": item.ic,
                "valid": item.valid,
                "sources": baseline_sources[item.ranks],
                "converged": item.fit.converged,
                "stationarity_residual": item.fit.stationarity_residual,
                "max_envelope_ratio": item.fit.max_envelope_ratio,
                "invalid_reasons": item.invalid_reasons,
                **item.start_diagnostics,
            }
            for item in all_records
            if item.ranks in baseline_candidate_vectors
        ],
        "cap_pilot_converged": cap_fit.converged,
        "cap_pilot_stationarity_residual": cap_fit.stationarity_residual,
        "cap_pilot_max_envelope_ratio": cap_fit.max_envelope_ratio,
        "cap_pilot": cap_pilot_diagnostics,
        "sensitivities": sensitivity,
    }
    return RankSelectionResult(selected, evaluated, preliminary, cap_fit, diagnostics)
