import numpy as np
import pytest

from dynamic_panel_econ.core import Coefficients, add
from dynamic_panel_econ.targets import (
    embedded_restriction,
    paper_index,
    target_direction,
    target_value,
)


def template():
    a = np.arange(48, dtype=float).reshape(8, 6)
    b = 2 * a
    return Coefficients([a], [b], np.zeros_like(a))


def test_paper_index_conversion():
    assert paper_index(8, 4) == 1
    assert paper_index(6, 2) == 2


def test_group_truth_generic_inner_product_matches_direct_formula():
    theta = template()
    groups = np.repeat([0, 1], 4)
    target = target_direction("A_G2_minus_G1_fixed_time", theta, groups)
    t0 = paper_index(6, 2)
    direct = theta.A[0][groups == 1, t0].mean() - theta.A[0][groups == 0, t0].mean()
    assert target_value(target.direction, theta) == pytest.approx(direct)


def test_restricted_targets_reassemble_exactly():
    theta = template()
    direction = target_direction("A_full_mean", theta).direction
    rows = np.arange(8)
    first = embedded_restriction(direction, rows, np.array([0, 2, 4]))
    second = embedded_restriction(direction, rows, np.array([1, 3, 5]))
    rebuilt = add(first, second)
    for actual, expected in zip(rebuilt.matrices(), direction.matrices(), strict=True):
        np.testing.assert_array_equal(actual, expected)
