from dynamic_panel_econ.rank_selection import model_dimension, one_coordinate_neighbors


def test_dimension_formula():
    assert model_dimension((1, 2, 0), 10, 8) == 1 * 17 + 2 * 16


def test_one_coordinate_candidates_are_valid_and_unique():
    neighbors = one_coordinate_neighbors((1, 0, 2), (3, 3, 3))
    assert neighbors == {(0, 0, 2), (2, 0, 2), (1, 1, 2), (1, 0, 1), (1, 0, 3)}
    assert len(neighbors) == len(set(neighbors))
