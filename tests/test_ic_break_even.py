from math import inf

import pytest

from dynamic_panel_econ.ic_break_even import (
    ICPoint,
    dimension_penalty,
    optimality_interval,
    pairwise_break_even,
    rank_increment_dimension,
    select_ic_candidate,
)


def point(ranks, qhat, dimension):
    return ICPoint(ranks, qhat, dimension)


def test_pairwise_break_even_and_inequality_direction() -> None:
    truth = point((1, 1, 1), 1.0, 30)
    underfit = point((0, 0, 0), 2.0, 0)
    c_star = pairwise_break_even(truth, underfit, base_penalty_per_dimension=2.0)
    assert c_star == pytest.approx(__import__("math").log(2.0) / 60)
    base = 2.0
    assert select_ic_candidate([truth, underfit], c_star / 2, base_penalty_per_dimension=base) == truth
    assert select_ic_candidate([truth, underfit], c_star * 2, base_penalty_per_dimension=base) == underfit


def test_lower_and_upper_bounds_intersect() -> None:
    truth = point((1, 1, 1), 1.0, 30)
    overfit = point((2, 1, 1), 0.8, 40)
    underfit = point((0, 1, 1), 1.5, 20)
    interval = optimality_interval(truth, [truth, overfit, underfit], base_penalty_per_dimension=1.0)
    assert interval.lower == pytest.approx(-__import__("math").log(0.8) / 10)
    assert interval.upper == pytest.approx(__import__("math").log(1.5) / 10)
    assert not interval.empty


def test_infinite_upper_bound() -> None:
    truth = point((1, 1, 1), 1.0, 30)
    overfit = point((2, 1, 1), 0.8, 40)
    interval = optimality_interval(truth, [overfit], base_penalty_per_dimension=1.0)
    assert interval.upper == inf
    assert not interval.empty


def test_empty_interval() -> None:
    truth = point((1, 1, 1), 1.0, 30)
    overfit = point((2, 1, 1), 0.5, 40)
    underfit = point((0, 1, 1), 1.1, 20)
    interval = optimality_interval(truth, [overfit, underfit], base_penalty_per_dimension=1.0)
    assert interval.empty
    assert interval.lower > interval.upper


def test_equal_dimension_better_competitor_makes_interval_empty() -> None:
    truth = point((1, 1, 1), 2.0, 30)
    competitor = point((0, 2, 1), 1.0, 30)
    assert optimality_interval(truth, [competitor], base_penalty_per_dimension=1.0).empty


def test_grid_selection_uses_dimension_and_rank_ties() -> None:
    first = point((1, 0, 0), 1.0, 10)
    second = point((0, 1, 0), 1.0, 10)
    assert select_ic_candidate([first, second], 1e-5, base_penalty_per_dimension=2.0) == second


def test_dimension_penalty_matches_rank_formula() -> None:
    assert rank_increment_dimension(0, 50, 50) == 99
    assert rank_increment_dimension(1, 50, 50) == 97
    assert dimension_penalty((1, 1, 1), 50, 50) > dimension_penalty((1, 0, 1), 50, 50)


def test_invalid_inputs_are_rejected() -> None:
    truth = point((1, 1, 1), 1.0, 30)
    with pytest.raises(ValueError):
        optimality_interval(truth, [], base_penalty_per_dimension=0.0)
    with pytest.raises(ValueError):
        select_ic_candidate([], 1.0, base_penalty_per_dimension=1.0)
    with pytest.raises(ValueError):
        rank_increment_dimension(-1, 50, 50)
