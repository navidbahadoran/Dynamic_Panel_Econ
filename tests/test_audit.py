from __future__ import annotations

import json
import subprocess
import sys
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import dynamic_panel_econ.rank_selection as rank_module
from dynamic_panel_econ.config import load_config
from dynamic_panel_econ.core import Coefficients, Design
from dynamic_panel_econ.dgp import DGP4_TRUTH_NAMES
from dynamic_panel_econ.estimation import FactorBlock, FitResult, NuclearFit
from dynamic_panel_econ.monte_carlo import (
    calibrate_design,
    classify_inference_status,
    run_monte_carlo,
    run_replication,
)
from dynamic_panel_econ.rank_selection import revision8_kappa, select_ranks
from dynamic_panel_econ.reporting import aggregate_run, make_tables


def _theta(ranks: tuple[int, ...], n: int = 6, t: int = 6) -> Coefficients:
    rng = np.random.default_rng(sum(ranks) + 10)
    matrices = []
    for rank in ranks:
        if rank == 0:
            matrices.append(np.zeros((n, t)))
        else:
            matrices.append(rng.normal(size=(n, rank)) @ rng.normal(size=(rank, t)))
    return Coefficients([matrices[0]], [matrices[1]], matrices[2])


def _fit(ranks: tuple[int, ...], objective: float = 1.0) -> FitResult:
    theta = _theta(ranks)
    blocks = [FactorBlock(np.empty((6, 0)), np.empty((6, 0))) for _ in ranks]
    return FitResult(theta, ranks, objective, True, 1, [objective], 0.0, 0.01, blocks)


def test_fresh_process_imports_installed_package():
    result = subprocess.run(
        [sys.executable, "-c", "import dynamic_panel_econ; print(dynamic_panel_econ.__version__)"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip() == "0.1.0"


def test_revision8_ic_penalty_formula_exactly():
    n, t, eta, spatial_dimension = 50, 70, 4.0, 1
    nt = n * t
    b_nt = nt ** (1 / (8 + eta)) * np.log(nt)
    expected = b_nt**2 * np.log(nt) ** (spatial_dimension + 3)
    assert revision8_kappa(n, t, eta_for_penalty=eta, spatial_dimension=spatial_dimension) == pytest.approx(expected)


def test_invalid_low_objective_candidate_cannot_win(monkeypatch):
    design = Design([np.ones((6, 6))], [np.ones((6, 6))])
    y = np.ones((6, 6))
    cap_fit = _fit((1, 1, 1), 2.0)
    zero = Coefficients([np.zeros((6, 6))], [np.zeros((6, 6))], np.zeros((6, 6)))
    nuclear = NuclearFit(zero, 1.0, 1.0, True, 1, [1.0], [[], [], []])
    monkeypatch.setattr(rank_module, "nuclear_path", lambda *args, **kwargs: [nuclear])
    monkeypatch.setattr(rank_module, "fit_fixed_rank", lambda *args, **kwargs: cap_fit)
    monkeypatch.setattr(
        rank_module,
        "fit_rank_adaptive_cap_pilot",
        lambda *args, **kwargs: (
            cap_fit,
            {
                "numerical_rank_before_thresholding": (1, 1, 1),
                "objective_stability_pass": True,
            },
        ),
    )

    def fake_candidate(*args, **kwargs):
        ranks = args[2]
        diagnostics = {
            "start_objectives": [1.0, 1.0],
            "start_stationarity_residuals": [0.0, 0.0],
            "best_two_objective_gap": 0.0,
            "objective_stability_pass": True,
            "third_start_used": False,
        }
        if ranks == (0, 0, 0):
            return _fit(ranks, 0.0), False, ["not_converged"], diagnostics
        return _fit(ranks, 1.0 + 0.1 * sum(ranks)), False, [], diagnostics

    monkeypatch.setattr(rank_module, "_fit_candidate", fake_candidate)
    result = select_ranks(
        y,
        design,
        (1, 1, 1),
        nuclear_epsilon=0.9,
        compute_rank_sensitivities=False,
        coefficient_bound=100.0,
        stationarity_tol=1.0,
    )
    assert result.candidates[(0, 0, 0)].valid is False
    assert result.selected.ranks != (0, 0, 0)


def test_candidate_requires_a_stable_pair_of_valid_start_objectives(monkeypatch):
    objectives = iter([1.0, 1.01, 1.02])
    monkeypatch.setattr(
        rank_module,
        "fit_fixed_rank",
        lambda *args, **kwargs: _fit((1, 1, 1), next(objectives)),
    )
    design = Design([np.ones((6, 6))], [np.ones((6, 6))])
    fit, third, reasons, diagnostics = rank_module._fit_candidate(
        np.ones((6, 6)),
        design,
        (1, 1, 1),
        (_theta((1, 1, 1)), _theta((1, 1, 1))),
        11,
        {},
        1e-6,
        1.0,
    )
    assert fit.objective == 1.0
    assert third is True
    assert "objective_stability_failed" in reasons
    assert diagnostics["objective_stability_pass"] is False
    assert diagnostics["start_objectives"] == [1.0, 1.01, 1.02]


def test_cap_pilot_is_rank_at_most_cap_and_may_return_lower_rank():
    zero = Coefficients([np.zeros((6, 6))], [np.zeros((6, 6))], np.zeros((6, 6)))
    preliminary = [
        NuclearFit(zero, penalty, 0.0, True, 1, [0.0], [[], [], []])
        for penalty in (1.0, 0.5)
    ]
    fit, diagnostics = rank_module.fit_rank_adaptive_cap_pilot(
        np.zeros((6, 6)),
        Design([np.ones((6, 6))], [np.ones((6, 6))]),
        (2, 2, 2),
        preliminary,
        0.1,
        seed=1,
        fit_options={"coefficient_bound": 9.0},
        stationarity_tolerance=1e-6,
        start_objective_stability_tol=1e-6,
        improvement_tolerance=1e-7,
        removal_tolerance=1e-7,
        max_steps=2,
    )
    assert fit.ranks == (0, 0, 0)
    assert fit.ranks != (2, 2, 2)
    assert diagnostics["objective_stability_pass"] is True
    assert diagnostics["attempted_route_count"] == 4
    assert diagnostics["valid_route_count"] >= 2
    assert any(
        any(attempt["start_rank_vector"])
        and tuple(attempt["final_rank_vector"]) == (0, 0, 0)
        and attempt["final_valid"]
        for attempt in diagnostics["outer_start_attempts"]
    )


def test_all_rank_sensitivities_use_completion_and_adaptive_cap_pilot(monkeypatch):
    zero = Coefficients([np.zeros((6, 6))], [np.zeros((6, 6))], np.zeros((6, 6)))
    nuclear = NuclearFit(zero, 1.0, 0.0, True, 1, [0.0], [[], [], []])
    monkeypatch.setattr(rank_module, "nuclear_path", lambda *args, **kwargs: [nuclear])
    cap_calls = []

    def fake_cap(y, design, caps, preliminary, threshold, **kwargs):
        cap_calls.append(tuple(caps))
        return _fit((0, 0, 0), 1.0), {
            "algorithm": "rank_adaptive_at_most_cap",
            "numerical_rank_before_thresholding": (0, 0, 0),
            "objective_stability_pass": True,
        }

    def fake_candidate(y, design, ranks, starts, *args, **kwargs):
        objective = 1.0 + 0.01 * sum(ranks)
        diagnostics = {
            "start_objectives": [objective, objective],
            "start_stationarity_residuals": [0.0, 0.0],
            "start_valid": [True, True],
            "best_two_objective_gap": 0.0,
            "objective_stability_pass": True,
            "third_start_used": False,
        }
        return _fit(ranks, objective), False, [], diagnostics

    monkeypatch.setattr(rank_module, "fit_rank_adaptive_cap_pilot", fake_cap)
    monkeypatch.setattr(rank_module, "_fit_candidate", fake_candidate)
    result = select_ranks(
        np.zeros((6, 6)),
        Design([np.ones((6, 6))], [np.ones((6, 6))]),
        (1, 1, 1),
        coefficient_bound=9.0,
        stationarity_tol=1.0,
        threshold_sensitivity_multipliers=[1.0],
        ic_sensitivity_multipliers=[1.0],
        larger_rank_caps=[2, 2, 2],
        compute_rank_sensitivities=True,
        compute_dense_grid_sensitivity=True,
        compute_larger_cap_sensitivity=True,
    )
    sensitivity = result.diagnostics["sensitivities"]
    assert cap_calls == [(1, 1, 1), (1, 1, 1), (1, 1, 1), (2, 2, 2)]
    assert sensitivity["threshold_multipliers"]["1.0"]["local_completion_applied"]
    assert sensitivity["ic_multipliers"]["1.0"]["local_completion_applied"]
    assert sensitivity["dense_grid_local_completion_applied"]
    assert sensitivity["larger_cap_local_completion_applied"]
    assert sensitivity["larger_cap_pilot_algorithm"] == "rank_adaptive_at_most_cap"


def _split_result(overrides: dict | None = None, count: int = 4):
    base = {
        "converged": True,
        "stationarity_residual": 0.0,
        "max_envelope_ratio": 0.1,
        "rank_supported": True,
        "target_supported": True,
        "riesz_converged": True,
        "riesz_target_stable": True,
    }
    splits = [dict(base) for _ in range(count)]
    if overrides:
        splits[0].update(overrides)
    return SimpleNamespace(
        riesz=SimpleNamespace(converged=True, target_rayleigh_quotient=1.0),
        variance=1.0,
        corrected=True,
        diagnostics={"split_fits": splits},
    )


@pytest.mark.parametrize(
    ("overrides", "count", "expected"),
    [
        (None, 3, "split_fit_not_converged"),
        ({"rank_supported": False}, 4, "split_rank_loss"),
        ({"target_supported": False}, 4, "split_target_support_loss"),
        ({"riesz_target_stable": False}, 4, "split_riesz_target_instability"),
        (None, 4, "success"),
    ],
)
def test_split_failure_diagnostics(overrides, count, expected):
    config = load_config("configs/mc/smoke.toml")
    assert classify_inference_status(_split_result(overrides, count), config) == expected


def test_dgp4_report_fields_include_a_b_raw_and_scaled_truths():
    required = {
        "A_G1_fixed_time_true",
        "A_G2_minus_G1_time_average_true",
        "B_G1_fixed_time_true",
        "B_G2_minus_G1_fixed_time_true",
        "B_G1_time_average_true",
        "B_G2_minus_G1_time_average_true",
        "A_G2_minus_G1_fixed_time_raw_true",
    }
    assert required.issubset(DGP4_TRUTH_NAMES)


def test_rank_stress_runner_dispatches_configured_true_rank(monkeypatch):
    config = load_config("configs/mc/rank_stress_smoke.toml")
    seen = []

    def fake_generate(*args, **kwargs):
        seen.append(tuple(args[3]))
        raise RuntimeError("dispatch observed")

    monkeypatch.setattr("dynamic_panel_econ.monte_carlo.generate_rank_stress_panel", fake_generate)
    calibration = {"c_h": 1.0, "c_xi": 1.0}
    for rank in ((1, 1, 1), (2, 1, 1), (1, 0, 2)):
        rows = run_replication((1, 20, 20, 0, rank), config, calibration)
        assert rows[0]["status"] == "unexpected_exception"
        assert json.loads(rows[0]["true_rank_vector"]) == list(rank)
    assert seen == [(1, 1, 1), (2, 1, 1), (1, 0, 2)]


def test_rank_at_cap_stops_primary_inference(monkeypatch):
    config = load_config("configs/mc/smoke.toml")
    config["inference"]["targets"] = []
    start_diagnostics = {
        "start_objectives": [1.0, 1.0],
        "start_stationarity_residuals": [0.0, 0.0],
        "best_two_objective_gap": 0.0,
        "objective_stability_pass": True,
        "third_start_used": False,
    }
    selected = SimpleNamespace(
        ranks=(3, 1, 1),
        fit=_fit((3, 1, 1)),
        ic=1.0,
        start_diagnostics=start_diagnostics,
    )
    diagnostics = {
        "candidate_rank_vectors": [(3, 1, 1)],
        "baseline_candidate_records": [
            {"ranks": (3, 1, 1), "ic": 1.0, "sources": ["rank_cap_pilot"]}
        ],
        "rank_cap_thresholded_vector": (3, 1, 1),
        "candidate_count_initial": 1,
        "candidate_count_final": 1,
        "smallest_ic": 1.0,
        "second_smallest_ic": None,
        "ic_gap": None,
        "selected_rank_at_cap": True,
        "cap_pilot_converged": True,
        "cap_pilot_stationarity_residual": 0.0,
        "cap_pilot_max_envelope_ratio": 0.1,
        "cap_pilot": {
            "numerical_rank_before_thresholding": (3, 1, 1),
            "outer_start_rank_vectors": [(3, 1, 1), (3, 1, 1)],
            "outer_final_rank_vectors": [(3, 1, 1), (3, 1, 1)],
            "outer_objectives": [1.0, 1.0],
            "best_two_objective_gap": 0.0,
            "objective_stability_pass": True,
            "outer_start_attempts": [],
        },
        "sensitivities": {},
    }
    fake = SimpleNamespace(selected=selected, diagnostics=diagnostics)
    monkeypatch.setattr("dynamic_panel_econ.monte_carlo.select_ranks", lambda *a, **k: fake)
    rows = run_replication((1, 20, 20, 0, None), config, {"c_h": 1.0, "c_xi": 1.0})
    assert [row["record_type"] for row in rows] == ["rank", "failure"]
    assert rows[1]["status"] == "rank_at_cap"


def _parquet_tree(root: Path, subdir: str):
    import pandas as pd

    return pd.concat(
        (pd.read_parquet(path) for path in sorted((root / subdir).glob("*.parquet"))),
        ignore_index=True,
    )


def test_serial_two_worker_end_to_end_equality_and_broad_schema(tmp_path):
    config = load_config("configs/mc/smoke.toml")
    config["run"].update(
        {
            "name": "serial_equality",
            "dgps": [1],
            "cells": [[20, 20]],
            "replications": 2,
            "chunk_size": 2,
            "n_jobs": 1,
            "output_root": str(tmp_path),
        }
    )
    config["dgp"]["group_gap_pilot_draws"] = 1
    config["inference"]["targets"] = ["A_full_mean", "B_full_mean"]
    _, serial_root = run_monte_carlo(config, overwrite=True)
    parallel = deepcopy(config)
    parallel["run"]["name"] = "parallel_equality"
    parallel["run"]["n_jobs"] = 2
    _, parallel_root = run_monte_carlo(parallel, overwrite=True)

    serial_raw = _parquet_tree(serial_root, "raw").sort_values(["replication", "target"])
    parallel_raw = _parquet_tree(parallel_root, "raw").sort_values(["replication", "target"])
    stable_columns = [
        "dgp", "N", "T", "replication", "target", "status", "true_rank_vector",
        "truth", "plugin_estimate", "corrected_estimate", "plugin_variance",
        "corrected_variance", "phi_full", "phi_time_sum", "phi_unit_sum",
        "split_fit_count",
    ]
    for column in stable_columns:
        left, right = serial_raw[column].reset_index(drop=True), parallel_raw[column].reset_index(drop=True)
        if left.dtype.kind in "fc":
            np.testing.assert_allclose(left, right, equal_nan=True)
        else:
            assert left.tolist() == right.tolist()
    assert set(
        [
            "plugin_estimate", "corrected_estimate", "plugin_standard_error",
            "corrected_standard_error", "plugin_error", "corrected_error",
            "phi_full", "phi_time_sum", "phi_unit_sum",
        ]
    ).issubset(serial_raw.columns)
    assert (serial_raw["split_fit_count"] == 4).all()
    assert (serial_raw["split_coefficient_fit_count"] == 4).all()
    fit_diagnostics = pd.read_parquet(serial_root / "fit_diagnostics.parquet")
    split_fits = fit_diagnostics.loc[
        fit_diagnostics["fit_type"].str.contains("_split_", na=False)
    ]
    assert len(split_fits) == 8
    assert (pd.to_numeric(split_fits["runtime_seconds"]) > 0).all()
    assert split_fits["diagnostic_context"].str.contains("_split_").all()
    assert serial_raw["true_target_projection_ratio"].notna().all()
    assert set(serial_raw["target_applicability"]) == {"theorem_covered"}
    for encoded in serial_raw["split_diagnostics_json"]:
        for split in json.loads(encoded):
            positive = [
                item
                for item in split["split_rank_singular_values"]
                if item["supplied_rank"] > 0
            ]
            assert positive
            assert all("sigma_r" in item and "sigma_r_over_sigma_1" in item for item in positive)


def test_rank_stress_design_normalizes_zero_slope_scale():
    config = load_config("configs/mc/rank_stress_smoke.toml")
    calibration = calibrate_design(config)[(1, 20, 20, (1, 0, 2))]
    assert calibration.c_xi == 1.0
    assert calibration.target_r2 is None
    assert calibration.diagnostics["r2_scale_identified"] is False


def test_rank_stress_end_to_end_runs_feasible_true_vectors(tmp_path):
    config = load_config("configs/mc/rank_stress_smoke.toml")
    config["run"]["output_root"] = str(tmp_path)
    _, root = run_monte_carlo(config, overwrite=True)
    ranks = _parquet_tree(root, "rank")
    raw = _parquet_tree(root, "raw")
    observed = set(ranks.get("true_rank_vector", [])) | set(
        raw.get("true_rank_vector", [])
    )
    assert observed == {
        "[1, 1, 1]",
        "[2, 1, 1]",
        "[1, 0, 2]",
    }
    aggregate_run(root)
    table_paths = make_tables(root)
    assert root / "tables" / "tab_mc_rank.tex" in table_paths
    rank_tex = (root / "tables" / "tab_mc_rank.tex").read_text(encoding="utf-8")
    assert "\\begin{longtable}" in rank_tex
