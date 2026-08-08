"""Revision-8 rank screening, valid post-refits, IC selection, and stability checks."""

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


class RankSelectionFailure(RuntimeError):
    """Raised when no valid post-refit can determine the IC minimum."""


class RankPilotFailure(RankSelectionFailure):
    """Raised when the imposed rank-cap pilot fails a paper diagnostic."""

    def __init__(self, message: str, diagnostics: dict[str, Any] | None = None) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics or {}


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
) -> list[str]:
    reasons = []
    if not fit.converged:
        reasons.append("not_converged")
    if not np.isfinite(fit.stationarity_residual) or fit.stationarity_residual > stationarity_tolerance:
        reasons.append("stationarity_high")
    if not np.isfinite(fit.max_envelope_ratio) or fit.max_envelope_ratio >= 1.0:
        reasons.append("coefficient_bound_active")
    actual = _numerical_rank_vector(fit.theta)
    if actual != required_ranks:
        reasons.append(f"numerical_rank_support:{actual}")
    if not np.isfinite(fit.objective):
        reasons.append("nonfinite_objective")
    return reasons


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
) -> tuple[FitResult, bool, list[str], dict[str, Any]]:
    prepared_starts = [adapt_initial(start, ranks) for start in starts]
    start_preparation: list[dict[str, float | bool]] = []
    if start_envelope_fraction is not None:
        prepared = []
        for start in prepared_starts:
            rescaled, diagnostics = _rescale_cap_start(
                start,
                float(fit_options.get("coefficient_bound", 9.0)),
                start_envelope_fraction,
            )
            prepared.append(rescaled)
            start_preparation.append(diagnostics)
        prepared_starts = prepared
    first = fit_fixed_rank(
        y, design, ranks, initial=prepared_starts[0], seed=seed, **fit_options
    )
    second = fit_fixed_rank(
        y, design, ranks, initial=prepared_starts[1], seed=seed, **fit_options
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
        options.append(fit_fixed_rank(y, design, ranks, seed=seed, **fit_options))
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
    coefficient_bound = float(fit_options.get("coefficient_bound", 9.0))
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
            start_rank, (start_theta, second_theta)
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
                collapsed_rank, (fit.theta, start_theta)
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
                    ranks, (move_initial, closest)
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
    if len(ordered) < 2:
        failure_diagnostics = {
            "attempted_route_count": len(routes),
            "valid_route_count": len(ordered),
            "stable_route_count": 0,
            "objective_stability_pass": False,
            "outer_start_attempts": route_attempts,
        }
        raise RankPilotFailure(
            "rank-adaptive cap pilot has fewer than two valid outer starts: "
            f"{len(ordered)} of {len(routes)} routes valid",
            failure_diagnostics,
        )
    outer_gap = abs(ordered[0][0].fit.objective - ordered[1][0].fit.objective) / max(
        1.0, abs(ordered[0][0].fit.objective)
    )
    if outer_gap > start_objective_stability_tol:
        failure_diagnostics = {
            "attempted_route_count": len(routes),
            "valid_route_count": len(ordered),
            "stable_route_count": 1,
            "best_two_objective_gap": outer_gap,
            "objective_stability_pass": False,
            "outer_start_attempts": route_attempts,
        }
        raise RankPilotFailure(
            "rank-adaptive cap pilot objective stability failed: "
            f"best-two normalized gap={outer_gap:.6g}",
            failure_diagnostics,
        )
    chosen = ordered[0][0]
    best_objective = ordered[0][0].fit.objective
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
        "best_two_objective_gap": outer_gap,
        "objective_stability_pass": True,
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


def select_ranks(
    y: np.ndarray,
    design: Design,
    caps: RankVector,
    *,
    seed: int | np.random.SeedSequence = 0,
    coefficient_bound: float = 9.0,
    max_sweeps: int = 200,
    objective_rtol: float = 1e-8,
    stationarity_tol: float = 1e-6,
    lstsq_rcond: float = 1e-10,
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
    threshold_sensitivity_multipliers: list[float] | tuple[float, ...] = (0.5, 1.0, 2.0),
    ic_sensitivity_multipliers: list[float] | tuple[float, ...] = (0.5, 1.0, 2.0),
    larger_rank_caps: list[int] | tuple[int, ...] = (4, 4, 4),
    compute_rank_sensitivities: bool = True,
    compute_dense_grid_sensitivity: bool = True,
    compute_larger_cap_sensitivity: bool = True,
) -> RankSelectionResult:
    """Select one full-panel rank vector; invalid fits can never minimize the IC."""

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
