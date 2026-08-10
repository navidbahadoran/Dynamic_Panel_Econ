"""Audit scientific equivalence separately from computational scaling.

The scientific replay uses four locked Revision-9 realizations and is never
written to canonical Monte Carlo tables.  The scaling workload consists only
of deterministic fixed-rank numerical fits; it is not statistical evidence.
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing
import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_info, threadpool_limits

from dynamic_panel_econ.config import load_config, resolve_execution_workers
from dynamic_panel_econ.core import Coefficients, Design, fitted_values, max_abs
from dynamic_panel_econ.estimation import fit_fixed_rank
from dynamic_panel_econ.mc_accounting import semantic_replication_id
from dynamic_panel_econ.monte_carlo import (
    Task,
    _process_memory_metrics,
    calibrate_design,
    iter_outer_task_results,
)

SCIENCE_SETTINGS = (1, 4)
SCALING_SETTINGS = (1, 4, 8, 12, 14)
SCALING_TASK_COUNT = 20
SCALING_KERNEL_REPEATS = 12
SCIENCE_TASKS: tuple[Task, ...] = (
    (1, 50, 50, 2026080901, None),
    (2, 50, 50, 2026080901, None),
    (3, 100, 100, 2026080901, None),
    (4, 100, 100, 2026080901, None),
)
FLOAT_RTOL = 1e-10
FLOAT_ATOL = 1e-12
_SCALING_THREAD_CONTROLLER: Any = None


def _semantic_id(task: Task | tuple[Any, ...]) -> str:
    dgp, n, t, replication, stress_rank = task
    return semantic_replication_id(
        int(dgp), int(n), int(t), int(replication), tuple(stress_rank or (1, 1, 1))
    )


def _scientific_records(rows_by_task: dict[Task, list[dict[str, Any]]]) -> dict[str, Any]:
    records: dict[str, Any] = {}
    for task in SCIENCE_TASKS:
        rows = rows_by_task[task]
        semantic_id = _semantic_id(task)
        hashes = {
            str(row.get("dgp_realization_hash"))
            for row in rows
            if row.get("dgp_realization_hash") is not None
        }
        if len(hashes) != 1:
            raise AssertionError(f"expected one DGP hash for {semantic_id}: {hashes}")
        records[f"{semantic_id}|dgp"] = {"dgp_realization_hash": next(iter(hashes))}
        rank_rows = [row for row in rows if row["record_type"] == "rank"]
        if len(rank_rows) != 1:
            raise AssertionError(f"expected one rank row for {semantic_id}")
        rank_row = rank_rows[0]
        records[f"{semantic_id}|rank_selection"] = {
            "selected_rank_vector": rank_row.get("selected_rank_vector"),
            "candidate_coverage": rank_row.get("candidate_coverage"),
            "primary_status": rank_row.get("primary_status"),
        }
        for candidate in json.loads(rank_row["rank_diagnostics_json"])["candidate_records"]:
            ranks = tuple(int(value) for value in candidate["ranks"])
            objective = float(candidate["objective"])
            records[f"{semantic_id}|candidate|{ranks}"] = {
                "rank_vector": ranks,
                "objective": objective,
                "Q_hat": 2.0 * objective,
                "dimension": int(candidate["dimension"]),
                "IC": float(candidate["ic"]),
                "valid": bool(candidate["valid"]),
                "invalid_reasons": tuple(candidate.get("invalid_reasons", [])),
                "sources": tuple(sorted(candidate.get("sources", []))),
            }
        for row in rows:
            if row["record_type"] not in {"target", "failure"}:
                continue
            target = str(row.get("target") or "__failure__")
            records[f"{semantic_id}|{row['record_type']}|{target}"] = {
                "target": target,
                "truth": row.get("truth"),
                "estimate": row.get("estimate"),
                "plugin_estimate": row.get("plugin_estimate"),
                "corrected_estimate": row.get("corrected_estimate"),
                "status": row.get("status"),
                "primary_status": row.get("primary_status"),
                "point_estimate_valid": row.get("point_estimate_valid"),
                "inference_valid": row.get("inference_valid"),
                "failure_detail": row.get("failure_detail"),
            }
    return json.loads(json.dumps(records, sort_keys=True))


def _compare_values(left: Any, right: Any) -> tuple[bool, float, float]:
    if left is None or right is None:
        return left is right, 0.0, 0.0
    if isinstance(left, (float, np.floating)) or isinstance(right, (float, np.floating)):
        left_float, right_float = float(left), float(right)
        if math.isnan(left_float) and math.isnan(right_float):
            return True, 0.0, 0.0
        if math.isinf(left_float) or math.isinf(right_float):
            return left_float == right_float, 0.0, 0.0
        absolute = abs(left_float - right_float)
        relative = absolute / max(1.0, abs(left_float))
        return bool(
            np.isclose(left_float, right_float, rtol=FLOAT_RTOL, atol=FLOAT_ATOL)
        ), absolute, relative
    return left == right, 0.0, 0.0


def _compare_records(
    baseline: dict[str, Any], comparison: dict[str, Any], requested: int, effective: int
) -> dict[str, Any]:
    structure_equal = baseline.keys() == comparison.keys()
    mismatch_count = 0
    maximum_absolute = 0.0
    maximum_relative = 0.0
    compared_fields = 0
    if structure_equal:
        for record_key in baseline:
            left, right = baseline[record_key], comparison[record_key]
            if left.keys() != right.keys():
                structure_equal = False
                mismatch_count += 1
                continue
            for field in left:
                compared_fields += 1
                equal, absolute, relative = _compare_values(left[field], right[field])
                mismatch_count += int(not equal)
                maximum_absolute = max(maximum_absolute, absolute)
                maximum_relative = max(maximum_relative, relative)
    return {
        "baseline_requested_n_jobs": 1,
        "comparison_requested_n_jobs": requested,
        "comparison_effective_n_jobs": effective,
        "scientific_record_count": len(comparison),
        "compared_field_count": compared_fields,
        "structure_equal": structure_equal,
        "mismatch_count": mismatch_count,
        "maximum_absolute_numeric_difference": maximum_absolute,
        "maximum_relative_numeric_difference": maximum_relative,
        "scientific_equivalence_pass": structure_equal and mismatch_count == 0,
        "floating_rtol": FLOAT_RTOL,
        "floating_atol": FLOAT_ATOL,
    }


def _initialize_scaling_worker() -> None:
    global _SCALING_THREAD_CONTROLLER
    _SCALING_THREAD_CONTROLLER = threadpool_limits(limits=1)


def _rank_one_matrix(rng: np.random.Generator, n: int, t: int, envelope: float) -> np.ndarray:
    left = rng.normal(size=n)
    right = rng.normal(size=t)
    matrix = np.outer(left, right)
    return envelope * matrix / np.max(np.abs(matrix))


def _scaling_worker(task_index: int) -> dict[str, Any]:
    """Run deterministic interior and box-constrained fixed-rank kernels."""

    cpu_started = time.process_time()
    wall_started = time.perf_counter()
    before = threadpool_info()
    rng = np.random.default_rng(904_000 + task_index)
    n = t = (50, 75, 100, 125)[task_index % 4]
    design = Design([rng.normal(size=(n, t))], [rng.normal(size=(n, t))])
    truth = Coefficients(
        [_rank_one_matrix(rng, n, t, 0.20)],
        [_rank_one_matrix(rng, n, t, 0.18)],
        _rank_one_matrix(rng, n, t, 0.16),
    )
    y = fitted_values(truth, design) + 0.01 * rng.normal(size=(n, t))
    constrained_y = np.full((n, t), 1.25 + 0.01 * (task_index % 3))
    interior_objective_sum = 0.0
    constrained_objective_sum = 0.0
    for repeat in range(SCALING_KERNEL_REPEATS):
        interior = fit_fixed_rank(
            y,
            design,
            (1, 1, 1),
            seed=904_000 + 100 * task_index + repeat,
            max_sweeps=30,
            coefficient_bound=2.0,
            diagnostic_context="parallel_scaling_interior",
        )
        constrained = fit_fixed_rank(
            constrained_y,
            design,
            (0, 0, 1),
            seed=914_000 + 100 * task_index + repeat,
            max_sweeps=20,
            coefficient_bound=0.5,
            constrained_subproblem_max_iterations=80,
            diagnostic_context="parallel_scaling_box_constrained",
        )
        if not constrained.diagnostics["constrained_fallback_used"]:
            raise AssertionError("scaling task failed to exercise constrained estimation")
        interior_objective_sum += float(interior.objective)
        constrained_objective_sum += float(constrained.objective)
    current_rss, peak_rss = _process_memory_metrics()
    after = threadpool_info()
    return {
        "task_index": task_index,
        "dimension": n,
        "kernel_repeats": SCALING_KERNEL_REPEATS,
        "interior_objective_sum": interior_objective_sum,
        "interior_stationarity": float(interior.stationarity_residual),
        "interior_max_abs": max_abs(interior.theta),
        "interior_converged": bool(interior.converged),
        "constrained_objective_sum": constrained_objective_sum,
        "constrained_stationarity": float(constrained.stationarity_residual),
        "constrained_max_abs": max_abs(constrained.theta),
        "constrained_status": constrained.diagnostics["constrained_solver_status"],
        "constrained_fallback_used": True,
        "worker_pid": os.getpid(),
        "process_cpu_seconds": time.process_time() - cpu_started,
        "worker_wall_seconds": time.perf_counter() - wall_started,
        "threadpool_info_before": before,
        "threadpool_info_after": after,
        "current_rss_bytes": current_rss,
        "peak_rss_bytes": peak_rss,
    }


def _scaling_records(results: list[dict[str, Any]]) -> dict[str, Any]:
    excluded = {
        "worker_pid",
        "process_cpu_seconds",
        "worker_wall_seconds",
        "threadpool_info_before",
        "threadpool_info_after",
        "current_rss_bytes",
        "peak_rss_bytes",
    }
    return {
        f"scaling_task_{row['task_index']:02d}": {
            key: value for key, value in row.items() if key not in excluded
        }
        for row in sorted(results, key=lambda item: int(item["task_index"]))
    }


def _run_scaling(requested: int) -> tuple[list[dict[str, Any]], float]:
    effective = min(requested, SCALING_TASK_COUNT)
    started = time.perf_counter()
    if effective == 1:
        _initialize_scaling_worker()
        results = [_scaling_worker(index) for index in range(SCALING_TASK_COUNT)]
    else:
        with ProcessPoolExecutor(
            max_workers=effective,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=_initialize_scaling_worker,
        ) as executor:
            results = list(executor.map(_scaling_worker, range(SCALING_TASK_COUNT)))
    return sorted(results, key=lambda item: int(item["task_index"])), time.perf_counter() - started


def _thread_rows(
    part: str, requested: int, effective: int, metrics: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for metric in metrics:
        task_id = metric.get("task_index")
        if task_id is None:
            task_id = _semantic_id(tuple(metric["task"]))
        for phase in ("before", "after"):
            for info in metric[f"threadpool_info_{phase}"]:
                rows.append(
                    {
                        "benchmark_part": part,
                        "requested_n_jobs": requested,
                        "effective_n_jobs": effective,
                        "task_id": task_id,
                        "worker_pid": metric["worker_pid"],
                        "phase": phase,
                        "user_api": info.get("user_api"),
                        "internal_api": info.get("internal_api"),
                        "prefix": info.get("prefix"),
                        "version": info.get("version"),
                        "num_threads": info.get("num_threads"),
                        "one_native_thread_policy_pass": int(info.get("num_threads", 0)) <= 1,
                    }
                )
    return rows


def _memory_rows(
    part: str, requested: int, effective: int, metrics: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    rows = []
    for pid in sorted({int(metric["worker_pid"]) for metric in metrics}):
        group = [metric for metric in metrics if int(metric["worker_pid"]) == pid]
        current = [metric["current_rss_bytes"] for metric in group if metric["current_rss_bytes"] is not None]
        peaks = [metric["peak_rss_bytes"] for metric in group if metric["peak_rss_bytes"] is not None]
        rows.append(
            {
                "benchmark_part": part,
                "requested_n_jobs": requested,
                "effective_n_jobs": effective,
                "worker_pid": pid,
                "tasks_completed_by_worker": len(group),
                "maximum_observed_current_rss_bytes": max(current) if current else np.nan,
                "process_lifetime_peak_rss_bytes": max(peaks) if peaks else np.nan,
                "memory_measurement": "Windows GetProcessMemoryInfo" if os.name == "nt" else "resource.getrusage",
            }
        )
    return rows


def _benchmark_row(
    part: str,
    requested: int,
    effective: int,
    task_count: int,
    wall: float,
    metrics: list[dict[str, Any]],
    deterministic_pass: bool,
) -> dict[str, Any]:
    logical_processors = int(os.cpu_count() or 1)
    cpu_seconds = sum(float(metric["process_cpu_seconds"]) for metric in metrics)
    per_pid_peaks = {
        int(metric["worker_pid"]): max(
            int(item["peak_rss_bytes"] or 0)
            for item in metrics
            if int(item["worker_pid"]) == int(metric["worker_pid"])
        )
        for metric in metrics
    }
    return {
        "benchmark_part": part,
        "requested_n_jobs": requested,
        "effective_n_jobs": effective,
        "number_of_tasks": task_count,
        "wall_clock_seconds": wall,
        "aggregate_process_cpu_seconds": cpu_seconds,
        "mean_total_cpu_utilization_percent": 100.0 * cpu_seconds / (wall * logical_processors),
        "peak_cpu_utilization_percent": np.nan,
        "peak_cpu_measurement_status": "not_measured_reliably",
        "peak_resident_memory_bytes_upper_bound": sum(per_pid_peaks.values()),
        "tasks_completed": task_count,
        "tasks_per_second": task_count / wall,
        "logical_processors": logical_processors,
        "deterministic_output_pass": deterministic_pass,
    }


def _report(
    benchmark: pd.DataFrame,
    equivalence: pd.DataFrame,
    threads: pd.DataFrame,
    memory: pd.DataFrame,
) -> str:
    science = benchmark[benchmark["benchmark_part"] == "scientific_equivalence"]
    scaling = benchmark[benchmark["benchmark_part"] == "parallel_scaling"]
    fastest = scaling.loc[scaling["wall_clock_seconds"].idxmin()]
    recommended = scaling[
        scaling["wall_clock_seconds"] <= 1.02 * fastest["wall_clock_seconds"]
    ].sort_values("requested_n_jobs").iloc[0]
    serial = scaling[scaling["requested_n_jobs"] == 1].iloc[0]
    speedup = float(serial["wall_clock_seconds"] / fastest["wall_clock_seconds"])
    thread_pass = bool(threads["one_native_thread_policy_pass"].all()) if len(threads) else False
    equivalence_pass = bool(equivalence["scientific_equivalence_pass"].all())
    memory_available = bool(memory["process_lifetime_peak_rss_bytes"].notna().any())
    lines = [
        "# Parallel execution audit",
        "",
        "The scientific section replays four locked Revision-9 realizations only. The scaling section is a deterministic computational workload using the same fixed-rank estimation kernel, including box-constrained estimation. Neither section is new Monte Carlo evidence.",
        "",
        "## Architecture",
        "",
        "Previously, pools were repeatedly created inside cell/chunk loops and only a small chunk could run concurrently. The revised executor uses one bounded Windows-spawn outer pool across DGP x N x T x semantic-replication tasks. Nuclear paths, cap pilots, candidate refits, local completion, and split fits remain sequential inside each worker; there is no nested process pool.",
        "",
        "Configuration and frozen calibrations are initialized once per worker. Panel and design arrays are generated inside their owning worker instead of copied from a parent cache.",
        "",
        "## A. Scientific equivalence",
        "",
        "| Requested | Effective | Tasks | Wall seconds | Tasks/s |",
        "|--:|--:|--:|--:|--:|",
    ]
    for row in science.itertuples():
        lines.append(f"| {row.requested_n_jobs} | {row.effective_n_jobs} | {row.number_of_tasks} | {row.wall_clock_seconds:.3f} | {row.tasks_per_second:.6g} |")
    lines += [
        "",
        f"Scientific equivalence passed: `{equivalence_pass}`. The comparison covers semantic IDs, realization hashes, candidate sets, selected ranks, objectives, Q_hat, IC values, target estimates, statuses, and failure classifications.",
        "",
        "A previously completed requested-8 replay also had effective worker count 4. Its checkpoint is preserved as `.superseded_scientific_njobs_8_effective4.json`, but its redundant timing is excluded from this redesigned benchmark and is not interpreted.",
        "",
        "## B. Parallel scaling",
        "",
        "| Requested | Effective | Tasks | Wall seconds | Tasks/s | Mean total CPU % | Peak RSS upper bound GiB |",
        "|--:|--:|--:|--:|--:|--:|--:|",
    ]
    for row in scaling.itertuples():
        lines.append(f"| {row.requested_n_jobs} | {row.effective_n_jobs} | {row.number_of_tasks} | {row.wall_clock_seconds:.3f} | {row.tasks_per_second:.6g} | {row.mean_total_cpu_utilization_percent:.2f} | {row.peak_resident_memory_bytes_upper_bound / 2**30:.3f} |")
    lines += [
        "",
        f"Fastest scaling setting: requested/effective `{int(fastest.requested_n_jobs)}`/`{int(fastest.effective_n_jobs)}`, with `{speedup:.2f}x` wall-time speedup over serial.",
        f"All scaling outputs deterministic across worker settings: `{bool(scaling['deterministic_output_pass'].all())}`.",
        f"One-native-thread policy passed for every detected BLAS/OpenMP library: `{thread_pass}`.",
        f"Per-process resident-memory measurement available: `{memory_available}`.",
        "",
        "Peak CPU was not sampled reliably; aggregate worker CPU time divided by wall time and 20 logical processors is reported as mean total CPU utilization. Peak RSS is a conservative sum of per-worker lifetime peaks, not a synchronized system-wide sample.",
        "",
        "## Recommendation",
        "",
        f"For this 14-physical-core/20-logical-processor/approximately-40-GB machine, use `--n-jobs {int(recommended.requested_n_jobs)}` for workloads with enough outer tasks. It is the smallest setting within 2% of the fastest wall time and therefore avoids extra worker memory for a negligible timing difference. Retain `--n-jobs 1` for deterministic debugging. Recheck Task Manager during the first separately authorized production launch because full selected-rank tasks are heavier than this scaling kernel.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("configs/mc/preflight_revision9_locked_selected.toml"))
    parser.add_argument("--output-root", type=Path, default=Path("results/mc/audit/parallel_execution"))
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    output_files = {
        "benchmark": args.output_root / "benchmark_results.csv",
        "equivalence": args.output_root / "scientific_equivalence.csv",
        "threads": args.output_root / "threadpool_diagnostics.csv",
        "memory": args.output_root / "memory_diagnostics.csv",
    }
    prior_benchmark = pd.read_csv(output_files["benchmark"]) if args.resume and output_files["benchmark"].exists() else pd.DataFrame()
    prior_equivalence = pd.read_csv(output_files["equivalence"]) if args.resume and output_files["equivalence"].exists() else pd.DataFrame()
    prior_threads = pd.read_csv(output_files["threads"]) if args.resume and output_files["threads"].exists() else pd.DataFrame()
    prior_memory = pd.read_csv(output_files["memory"]) if args.resume and output_files["memory"].exists() else pd.DataFrame()

    benchmark_rows: list[dict[str, Any]] = []
    equivalence_rows: list[dict[str, Any]] = []
    thread_rows: list[dict[str, Any]] = []
    memory_rows: list[dict[str, Any]] = []

    baseline_checkpoint = args.output_root / ".scientific_njobs_1.json"
    baseline_records = json.loads(baseline_checkpoint.read_text(encoding="utf-8")) if args.resume and baseline_checkpoint.exists() else None
    if not prior_benchmark.empty:
        for row in prior_benchmark.to_dict("records"):
            requested = int(row["requested_n_jobs"])
            if requested not in SCIENCE_SETTINGS:
                continue
            benchmark_rows.append({
                "benchmark_part": "scientific_equivalence",
                "requested_n_jobs": requested,
                "effective_n_jobs": int(row["effective_n_jobs"]),
                "number_of_tasks": int(row.get("outer_tasks_available", row.get("number_of_tasks", 4))),
                "wall_clock_seconds": float(row["wall_clock_seconds"]),
                "aggregate_process_cpu_seconds": float(row["aggregate_process_cpu_seconds"]),
                "mean_total_cpu_utilization_percent": float(row["mean_total_cpu_utilization_percent"]),
                "peak_cpu_utilization_percent": np.nan,
                "peak_cpu_measurement_status": "not_measured_reliably",
                "peak_resident_memory_bytes_upper_bound": float(row["peak_resident_memory_bytes_upper_bound"]),
                "tasks_completed": int(row["tasks_completed"]),
                "tasks_per_second": float(row.get("effective_tasks_per_second", row.get("tasks_per_second"))),
                "logical_processors": int(row["logical_processors"]),
                "deterministic_output_pass": True,
            })
    if not prior_equivalence.empty:
        equivalence_rows = [row for row in prior_equivalence.to_dict("records") if int(row["comparison_requested_n_jobs"]) in SCIENCE_SETTINGS]
    if not prior_threads.empty:
        thread_rows = [{"benchmark_part": "scientific_equivalence", **row} for row in prior_threads.to_dict("records") if int(row["requested_n_jobs"]) in SCIENCE_SETTINGS]
    if not prior_memory.empty:
        memory_rows = [{"benchmark_part": "scientific_equivalence", **row} for row in prior_memory.to_dict("records") if int(row["requested_n_jobs"]) in SCIENCE_SETTINGS]

    completed_science = {int(row["requested_n_jobs"]) for row in benchmark_rows}
    if completed_science != set(SCIENCE_SETTINGS):
        base_config = load_config(args.config)
        base_config["run"]["blas_threads"] = 1
        calibrations = {key: asdict(value) for key, value in calibrate_design(base_config).items()}
        required_keys = {(task[0], task[1], task[2], task[4]) for task in SCIENCE_TASKS}
        if not required_keys <= calibrations.keys():
            raise AssertionError("frozen calibration is missing a scientific-equivalence task")
        for requested in SCIENCE_SETTINGS:
            if requested in completed_science:
                continue
            config = json.loads(json.dumps(base_config))
            resolve_execution_workers(config, requested_n_jobs=requested)
            effective = min(requested, len(SCIENCE_TASKS))
            rows_by_task: dict[Task, list[dict[str, Any]]] = {}
            metrics: list[dict[str, Any]] = []
            started = time.perf_counter()
            for task, rows, metric in iter_outer_task_results(list(SCIENCE_TASKS), config, calibrations, effective_n_jobs=effective, collect_metrics=True):
                rows_by_task[task] = rows
                if metric is not None:
                    metrics.append(metric)
            wall = time.perf_counter() - started
            scientific = _scientific_records(rows_by_task)
            (args.output_root / f".scientific_njobs_{requested}.json").write_text(json.dumps(scientific, sort_keys=True, allow_nan=True), encoding="utf-8")
            if baseline_records is None:
                baseline_records = scientific
            equivalence_rows.append(_compare_records(baseline_records, scientific, requested, effective))
            benchmark_rows.append(_benchmark_row("scientific_equivalence", requested, effective, len(SCIENCE_TASKS), wall, metrics, True))
            thread_rows.extend(_thread_rows("scientific_equivalence", requested, effective, metrics))
            memory_rows.extend(_memory_rows("scientific_equivalence", requested, effective, metrics))

    scaling_baseline: dict[str, Any] | None = None
    for requested in SCALING_SETTINGS:
        results, wall = _run_scaling(requested)
        records = _scaling_records(results)
        if scaling_baseline is None:
            scaling_baseline = records
        comparison = _compare_records(scaling_baseline, records, requested, min(requested, SCALING_TASK_COUNT))
        if not comparison["scientific_equivalence_pass"]:
            raise AssertionError("parallel scaling changed deterministic kernel outputs")
        benchmark_rows.append(_benchmark_row("parallel_scaling", requested, min(requested, SCALING_TASK_COUNT), SCALING_TASK_COUNT, wall, results, True))
        thread_rows.extend(_thread_rows("parallel_scaling", requested, min(requested, SCALING_TASK_COUNT), results))
        memory_rows.extend(_memory_rows("parallel_scaling", requested, min(requested, SCALING_TASK_COUNT), results))
        print(f"completed scaling requested_n_jobs={requested}, effective_n_jobs={min(requested, SCALING_TASK_COUNT)}, wall={wall:.3f}s", flush=True)

    benchmark = pd.DataFrame(benchmark_rows)
    benchmark["_part_order"] = benchmark["benchmark_part"].map(
        {"scientific_equivalence": 0, "parallel_scaling": 1}
    )
    benchmark = (
        benchmark.sort_values(["_part_order", "requested_n_jobs"])
        .drop(columns="_part_order")
        .reset_index(drop=True)
    )
    equivalence = pd.DataFrame(equivalence_rows).sort_values("comparison_requested_n_jobs").reset_index(drop=True)
    threads = pd.DataFrame(thread_rows).sort_values(["benchmark_part", "requested_n_jobs", "worker_pid", "task_id"], na_position="last").reset_index(drop=True)
    memory = pd.DataFrame(memory_rows).sort_values(["benchmark_part", "requested_n_jobs", "worker_pid"]).reset_index(drop=True)
    if not equivalence["scientific_equivalence_pass"].all():
        raise AssertionError("parallel replay changed a scientific output")
    if len(threads) and not threads["one_native_thread_policy_pass"].all():
        raise AssertionError("a detected worker-native threadpool exceeded one thread")
    benchmark.to_csv(output_files["benchmark"], index=False)
    equivalence.to_csv(output_files["equivalence"], index=False)
    threads.to_csv(output_files["threads"], index=False)
    memory.to_csv(output_files["memory"], index=False)
    (args.output_root / "parallel_execution_report.md").write_text(_report(benchmark, equivalence, threads, memory), encoding="utf-8")


if __name__ == "__main__":
    main()
