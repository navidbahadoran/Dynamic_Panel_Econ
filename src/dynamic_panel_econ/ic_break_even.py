"""Offline break-even diagnostics for the Revision-8 rank-selection IC."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from math import inf, log

import numpy as np

from .rank_selection import model_dimension, revision8_kappa

RankVector = tuple[int, ...]


@dataclass(frozen=True, slots=True)
class ICPoint:
    """The objective and dimension of one valid, already-fitted rank candidate."""

    ranks: RankVector
    qhat: float
    dimension: int

    @property
    def log_qhat(self) -> float:
        return log(max(self.qhat, np.finfo(float).tiny))


@dataclass(frozen=True, slots=True)
class BreakEvenInterval:
    """Nonnegative multipliers for which the reference candidate is IC-optimal."""

    lower: float
    upper: float
    empty: bool
    lower_binding_rank: RankVector | None = None
    upper_binding_rank: RankVector | None = None


def dimension_penalty(
    ranks: RankVector,
    n: int,
    t: int,
    *,
    eta_for_penalty: float = 4.0,
    spatial_dimension: int = 1,
) -> float:
    """Return ``kappa_base * d(r) / (N*T)`` before applying ``c_kappa``."""

    return revision8_kappa(
        n,
        t,
        eta_for_penalty=eta_for_penalty,
        spatial_dimension=spatial_dimension,
    ) * model_dimension(ranks, n, t) / (n * t)


def pairwise_break_even(
    reference: ICPoint,
    competitor: ICPoint,
    *,
    base_penalty_per_dimension: float = 1.0,
) -> float:
    """Solve the equality in the affine IC comparison for ``c_kappa``.

    The reference wins when
    ``c * (P_reference - P_competitor) <= log(Q_competitor/Q_reference)``.
    """

    if base_penalty_per_dimension <= 0:
        raise ValueError("base_penalty_per_dimension must be positive")
    penalty_difference = base_penalty_per_dimension * (
        reference.dimension - competitor.dimension
    )
    if penalty_difference == 0:
        return inf if reference.log_qhat <= competitor.log_qhat else -inf
    return (competitor.log_qhat - reference.log_qhat) / penalty_difference


def optimality_interval(
    reference: ICPoint,
    competitors: Iterable[ICPoint],
    *,
    base_penalty_per_dimension: float,
) -> BreakEvenInterval:
    """Intersect all pairwise IC inequalities over ``c_kappa >= 0``."""

    if base_penalty_per_dimension <= 0:
        raise ValueError("base_penalty_per_dimension must be positive")
    lower = 0.0
    upper = inf
    lower_rank: RankVector | None = None
    upper_rank: RankVector | None = None
    for competitor in competitors:
        if competitor.ranks == reference.ranks:
            continue
        delta_penalty = base_penalty_per_dimension * (
            reference.dimension - competitor.dimension
        )
        delta_log_q = competitor.log_qhat - reference.log_qhat
        if delta_penalty > 0:
            bound = delta_log_q / delta_penalty
            if bound < upper:
                upper, upper_rank = bound, competitor.ranks
        elif delta_penalty < 0:
            bound = delta_log_q / delta_penalty
            if bound > lower:
                lower, lower_rank = bound, competitor.ranks
        elif delta_log_q < 0:
            return BreakEvenInterval(0.0, -inf, True, None, competitor.ranks)
    return BreakEvenInterval(
        lower=lower,
        upper=upper,
        empty=upper < max(0.0, lower),
        lower_binding_rank=lower_rank,
        upper_binding_rank=upper_rank,
    )


def select_ic_candidate(
    candidates: Sequence[ICPoint],
    multiplier: float,
    *,
    base_penalty_per_dimension: float,
) -> ICPoint:
    """Reproduce the rank selector's deterministic IC/dimension/rank tie-breaking."""

    if not candidates:
        raise ValueError("at least one valid candidate is required")
    if multiplier < 0:
        raise ValueError("multiplier must be nonnegative")
    return min(
        candidates,
        key=lambda point: (
            point.log_qhat + multiplier * base_penalty_per_dimension * point.dimension,
            point.dimension,
            point.ranks,
        ),
    )


def rank_increment_dimension(rank: int, n: int, t: int) -> int:
    """Dimension increment for one matrix when its rank rises from r to r+1."""

    if rank < 0:
        raise ValueError("rank must be nonnegative")
    return n + t - 2 * rank - 1
