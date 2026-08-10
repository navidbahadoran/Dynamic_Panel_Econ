from __future__ import annotations

import inspect
import pickle
from copy import deepcopy
from pathlib import Path

import pandas as pd

from dynamic_panel_econ import monte_carlo
from dynamic_panel_econ.cli import build_run_parser, resolve_run_args, resolved_config_text
from dynamic_panel_econ.config import load_config, resolve_execution_workers


def test_cli_n_jobs_and_precedence() -> None:
    parser = build_run_parser()
    args = parser.parse_args(
        ["--config", "configs/mc/pilot.toml", "--n-jobs", "12", "--dry-run"]
    )
    config = resolve_run_args(args)
    assert config["run"]["n_jobs"] == 12
    assert config["run"]["requested_n_jobs"] == 12
    assert config["run"]["effective_n_jobs"] == min(
        12, config["run"]["number_of_outer_tasks"]
    )
    assert "requested_n_jobs = 12" in resolved_config_text(config)
    assert f"effective_n_jobs = {config['run']['effective_n_jobs']}" in resolved_config_text(
        config
    )


def test_toml_n_jobs_and_effective_task_cap() -> None:
    config = load_config("configs/mc/pilot.toml")
    assert config["run"]["n_jobs"] == 4
    assert config["run"]["requested_n_jobs"] == 4
    tiny = deepcopy(config)
    tiny["run"].update({"dgps": [1], "cells": [[10, 10]], "replications": 1})
    resolve_execution_workers(tiny, requested_n_jobs=14)
    assert tiny["run"]["requested_n_jobs"] == 14
    assert tiny["run"]["effective_n_jobs"] == 1
    assert tiny["run"]["number_of_outer_tasks"] == 1


def test_n_jobs_one_does_not_create_process_pool(monkeypatch) -> None:
    task = (1, 10, 10, 0, None)
    config = load_config("configs/mc/smoke.toml")
    initialized: list[int] = []

    def fake_initialize(_config, _calibrations, native_threads):
        initialized.append(native_threads)

    monkeypatch.setattr(monte_carlo, "_initialize_outer_worker", fake_initialize)
    monkeypatch.setattr(monte_carlo, "_outer_worker", lambda saved: [{"task": saved}])
    monkeypatch.setattr(
        monte_carlo,
        "ProcessPoolExecutor",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("pool created")),
    )
    output = list(
        monte_carlo.iter_outer_task_results(
            [task], config, {}, effective_n_jobs=1
        )
    )
    assert initialized == [1]
    assert output == [(task, [{"task": task}], None)]


def test_worker_entry_point_is_pickle_safe_and_run_script_is_guarded() -> None:
    pickle.dumps(monte_carlo._outer_worker)
    source = Path("scripts/run_mc.py").read_text(encoding="utf-8")
    assert 'if __name__ == "__main__":' in source
    assert "multiprocessing.get_context(\"spawn\")" in inspect.getsource(
        monte_carlo.iter_outer_task_results
    )


def test_no_nested_process_pool_and_thread_limit_is_explicit(monkeypatch) -> None:
    source = inspect.getsource(monte_carlo)
    assert source.count("with ProcessPoolExecutor(") == 1
    calls: list[int] = []
    sentinel = object()

    def fake_limits(*, limits):
        calls.append(limits)
        return sentinel

    monkeypatch.setattr(monte_carlo, "threadpool_limits", fake_limits)
    monte_carlo._initialize_outer_worker({}, {}, 1)
    assert calls == [1]
    assert monte_carlo._WORKER_THREAD_LIMIT_CONTROLLER is sentinel


def test_deterministic_output_row_order() -> None:
    rows = [
        {"semantic_replication_id": "b", "record_type": "target", "target": "B"},
        {"semantic_replication_id": "a", "record_type": "target", "target": "B"},
        {"semantic_replication_id": "a", "record_type": "target", "target": "A"},
    ]
    ordered = sorted(rows, key=monte_carlo._deterministic_row_key)
    assert [(row["semantic_replication_id"], row["target"]) for row in ordered] == [
        ("a", "A"),
        ("a", "B"),
        ("b", "B"),
    ]


def test_locked_preflight_fixed_selected_hashes_remain_matched() -> None:
    path = Path("results/mc/preflight_revision9_locked/matched_dgp_hashes.csv")
    matched = pd.read_csv(path)
    assert len(matched) == 24
    assert matched["hashes_match"].all()
    assert (matched["fixed_rank"] == matched["selected_rank"]).all()
