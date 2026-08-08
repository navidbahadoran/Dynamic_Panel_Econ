from dynamic_panel_econ.config import config_hash, load_config
from dynamic_panel_econ.dgp import generate_panel
from dynamic_panel_econ.seeds import seed_sequence


def test_stable_seed_is_independent_of_call_order():
    first = generate_panel(2, 8, 8, seed_sequence(7, "dgp", 2, 8, 8, 3))
    generate_panel(1, 8, 8, seed_sequence(7, "unrelated"))
    second = generate_panel(2, 8, 8, seed_sequence(7, "dgp", 2, 8, 8, 3))
    assert (first.y == second.y).all()


def test_config_hash_changes_with_statistical_setting():
    config = load_config("configs/mc/smoke.toml")
    original = config_hash(config)
    config["dgp"]["rho_s"] = 0.4
    assert config_hash(config) != original
