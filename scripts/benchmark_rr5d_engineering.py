"""Run the RR5d non-scientific deterministic solver/parallel benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import os
import pickle
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_info, threadpool_limits

from dynamic_panel_econ.cap_plus_one import fit_cap_plus_one
from dynamic_panel_econ.core import Coefficients, Design, fitted_values
from dynamic_panel_econ.estimation import fit_fixed_rank
from dynamic_panel_econ.monte_carlo import _process_memory_metrics

WORKER_COUNTS = (1, 4, 8, 12, 14)
TASK_COUNT = 24
_THREAD_CONTROLLER: Any = None


def _initialize_worker() -> None:
    global _THREAD_CONTROLLER
    _THREAD_CONTROLLER = threadpool_limits(limits=1)


def _fixture(task_index: int) -> tuple[np.ndarray, Design, Coefficients, tuple[int, int]]:
    shapes = ((50, 50), (100, 100), (100, 200), (200, 100))
    n, t = shapes[task_index % len(shapes)]
    rng = np.random.default_rng(900_000 + task_index)
    lag = rng.normal(size=(n, t))
    covariate = rng.normal(size=(n, t))
    ranks = (1 + task_index % 2, 1, 1 + (task_index // 2) % 2)
    matrices = [
        rng.normal(size=(n, rank)) @ rng.normal(size=(rank, t)) / (8.0 * rank)
        for rank in ranks
    ]
    theta = Coefficients([matrices[0]], [matrices[1]], matrices[2])
    y = fitted_values(theta, Design([lag], [covariate]))
    if task_index % 6 == 0:
        y = y + 0.02 * np.outer(np.linspace(-1, 1, n), np.linspace(1, -1, t))
    return y, Design([lag], [covariate]), theta, (n, t)


def _scientific_hash(fit: Any) -> str:
    digest = hashlib.sha256()
    for matrix in fit.theta.matrices():
        array = np.ascontiguousarray(matrix, dtype=np.float64)
        digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
        digest.update(array.tobytes())
    digest.update(np.float64(fit.objective).tobytes())
    digest.update(np.float64(fit.stationarity_residual).tobytes())
    digest.update(str(bool(fit.converged)).encode("ascii"))
    return digest.hexdigest()


def _benchmark_task(task_index: int) -> dict[str, Any]:
    cpu_started = time.process_time()
    wall_started = time.perf_counter()
    y, design, theta, shape = _fixture(task_index)
    pilot_started = time.perf_counter()
    pilot = fit_cap_plus_one(
        y,
        design,
        (4, 4, 4),
        seed=70_000 + task_index,
        max_sweeps=12,
        stationarity_tol=1e-6,
        constrained_kkt_tolerance=1e-4,
    )
    pilot_seconds = time.perf_counter() - pilot_started
    final_started = time.perf_counter()
    final = fit_fixed_rank(
        y,
        design,
        (1 + task_index % 2, 1, 1 + (task_index // 2) % 2),
        initial=theta,
        seed=80_000 + task_index,
        max_sweeps=12,
    )
    final_fit_seconds = time.perf_counter() - final_started
    current_rss, peak_rss = _process_memory_metrics()
    thread_info = threadpool_info()
    result = {
        "task_id": f"engineering_{task_index:03d}",
        "task_index": task_index,
        "N": shape[0],
        "T": shape[1],
        "pilot_hash": _scientific_hash(pilot),
        "final_hash": _scientific_hash(final),
        "pilot_objective": pilot.objective,
        "pilot_stationarity": pilot.stationarity_residual,
        "pilot_converged": pilot.converged,
        "pilot_status": pilot.diagnostics.get("constrained_solver_status"),
        "final_objective": final.objective,
        "final_stationarity": final.stationarity_residual,
        "final_converged": final.converged,
        "pilot_seconds": pilot_seconds,
        "final_fit_seconds": final_fit_seconds,
        "task_wall_seconds": time.perf_counter() - wall_started,
        "process_cpu_seconds": time.process_time() - cpu_started,
        "worker_pid": os.getpid(),
        "current_rss_bytes": current_rss,
        "peak_rss_bytes": peak_rss,
        "native_thread_counts": [
            int(item.get("num_threads", 0)) for item in thread_info
        ],
    }
    serialization_started = time.perf_counter()
    pickle.dumps(result, protocol=pickle.HIGHEST_PROTOCOL)
    result["serialization_seconds"] = time.perf_counter() - serialization_started
    return result


def _atomic_json(destination: Path, value: Any) -> float:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".partial")
    started = time.perf_counter()
    content = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    with temporary.open("wb") as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, destination)
    return time.perf_counter() - started


def _run_once(worker_count: int, repetition: int, output_root: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    tasks = list(range(TASK_COUNT))
    started = time.perf_counter()
    completions: list[float] = []
    results: list[dict[str, Any]] = []
    write_seconds = 0.0
    queue_depth: list[int] = []
    if worker_count == 1:
        _initialize_worker()
        for task in tasks:
            result = _benchmark_task(task)
            write_seconds += _atomic_json(
                output_root / f"w{worker_count}_rep{repetition}" / f"{result['task_id']}.json",
                result,
            )
            results.append(result)
            completions.append(time.perf_counter() - started)
            queue_depth.append(TASK_COUNT - len(results))
    else:
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=context,
            initializer=_initialize_worker,
        ) as executor:
            pending = {executor.submit(_benchmark_task, task): task for task in tasks}
            while pending:
                completed, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in completed:
                    pending.pop(future)
                    result = future.result()
                    write_seconds += _atomic_json(
                        output_root
                        / f"w{worker_count}_rep{repetition}"
                        / f"{result['task_id']}.json",
                        result,
                    )
                    results.append(result)
                    completions.append(time.perf_counter() - started)
                    queue_depth.append(len(pending))
    wall = time.perf_counter() - started
    results.sort(key=lambda item: item["task_id"])
    tail_index = max(0, len(completions) - worker_count)
    idle_tail = max(completions) - sorted(completions)[tail_index]
    logical_cpus = os.cpu_count() or 1
    cpu_seconds = sum(float(item["process_cpu_seconds"]) for item in results)
    task_seconds = sum(float(item["task_wall_seconds"]) for item in results)
    unique_peak = {}
    for item in results:
        pid = int(item["worker_pid"])
        unique_peak[pid] = max(unique_peak.get(pid, 0), int(item["peak_rss_bytes"] or 0))
    first_overhead = max(
        0.0,
        min(completions) - min(float(item["task_wall_seconds"]) for item in results),
    )
    summary = {
        "requested_n_jobs": worker_count,
        "effective_n_jobs": min(worker_count, TASK_COUNT),
        "repetition": repetition,
        "number_of_tasks": TASK_COUNT,
        "wall_seconds": wall,
        "tasks_per_minute": TASK_COUNT * 60.0 / wall,
        "process_cpu_seconds": cpu_seconds,
        "cpu_utilization_percent_of_logical_capacity": 100.0 * cpu_seconds / (wall * logical_cpus),
        "worker_utilization_percent": 100.0 * task_seconds / (wall * min(worker_count, TASK_COUNT)),
        "peak_worker_rss_bytes": max(unique_peak.values(), default=0),
        "summed_worker_peak_rss_bytes": sum(unique_peak.values()),
        "process_spawn_dispatch_overhead_seconds": first_overhead,
        "serialization_seconds": sum(float(item["serialization_seconds"]) for item in results),
        "output_write_seconds": write_seconds,
        "idle_tail_seconds": idle_tail,
        "maximum_queue_depth": max(queue_depth, default=0),
        "pilot_seconds": sum(float(item["pilot_seconds"]) for item in results),
        "final_fit_seconds": sum(float(item["final_fit_seconds"]) for item in results),
        "failures": 0,
        "pilot_nonconverged": sum(not bool(item["pilot_converged"]) for item in results),
        "max_native_threads": max(
            (max(item["native_thread_counts"], default=0) for item in results), default=0
        ),
    }
    return summary, results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("results/design/revision10_ridge_ratio/rr5d_engineering/benchmark_raw"),
    )
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    summaries: list[dict[str, Any]] = []
    baseline: dict[str, tuple[str, str]] | None = None
    equivalence_rows: list[dict[str, Any]] = []
    for repetition in range(args.repetitions):
        order = WORKER_COUNTS[repetition % len(WORKER_COUNTS) :] + WORKER_COUNTS[: repetition % len(WORKER_COUNTS)]
        for workers in order:
            summary, results = _run_once(workers, repetition + 1, args.output_root)
            signatures = {
                str(item["task_id"]): (str(item["pilot_hash"]), str(item["final_hash"]))
                for item in results
            }
            if workers == 1 and baseline is None:
                baseline = signatures
            equal = baseline is None or signatures == baseline
            equivalence_rows.append(
                {
                    "requested_n_jobs": workers,
                    "repetition": repetition + 1,
                    "scientific_outputs_equal_to_serial": equal,
                }
            )
            if not equal:
                raise AssertionError(f"worker-count scientific mismatch: {workers}")
            if int(summary["max_native_threads"]) > 1:
                raise AssertionError("native threadpool exceeded one thread")
            summaries.append(summary)
            print(
                f"workers={workers} repetition={repetition + 1} wall={summary['wall_seconds']:.3f}s",
                flush=True,
            )
    output = args.output_root.parent
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).sort_values(["requested_n_jobs", "repetition"]).to_csv(
        output / "deterministic_parallel_benchmark.csv", index=False, lineterminator="\n"
    )
    pd.DataFrame(equivalence_rows).sort_values(["requested_n_jobs", "repetition"]).to_csv(
        output / "parallel_equivalence.csv", index=False, lineterminator="\n"
    )


if __name__ == "__main__":
    main()
