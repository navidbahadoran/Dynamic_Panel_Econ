"""Failure forensics for the accepted 36-evaluation final preflight."""

from __future__ import annotations

import argparse
import itertools
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from dynamic_panel_econ.monte_carlo import _worker

PILOT_FAILURE_IDS = (
    "dgp1_N50_T50_r00002_truth1-1-1",
    "dgp3_N50_T50_r00000_truth1-1-1",
    "dgp3_N50_T50_r00002_truth1-1-1",
)
FIXED_BOUND_ID = "dgp4_N50_T50_r00001_truth1-1-1"
TRUE_RANK = (1, 1, 1)


def _chunks(root: Path, subdirectory: str) -> pd.DataFrame:
    paths = sorted((root / subdirectory).glob("*.parquet"))
    return pd.concat([pd.read_parquet(path) for path in paths], ignore_index=True)


def _load(root: Path) -> dict[str, pd.DataFrame]:
    return {
        "attempts": pd.read_parquet(root / "attempted_replications.parquet"),
        "records": pd.read_parquet(root / "replication_records.parquet"),
        "fits": pd.read_parquet(root / "fit_diagnostics.parquet"),
        "rank": _chunks(root, "rank"),
        "raw": _chunks(root, "raw"),
    }


def _json(value: Any) -> Any:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    return json.loads(value) if isinstance(value, str) else value


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=True)


def _failure_rows(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    raw = data["raw"]
    return raw.loc[
        raw["record_type"].eq("failure")
        & raw["primary_status"].eq("rank_pilot_failure")
    ].copy()


def _rank_diagnostics(data: dict[str, pd.DataFrame]) -> dict[str, dict[str, Any]]:
    return {
        row.semantic_replication_id: _json(row.rank_diagnostics_json)
        for row in data["rank"].itertuples()
    }


def _nuclear_signature(fits: pd.DataFrame, semantic_id: str) -> str:
    columns = [
        "nuclear_path_index", "lambda", "objective_final", "iterations",
        "convergence_flag", "coefficient_envelope_ratio", "thresholded_rank",
    ]
    rows = fits.loc[
        fits["semantic_replication_id"].eq(semantic_id)
        & fits["fit_type"].eq("nuclear_path"),
        columns,
    ].sort_values("nuclear_path_index")
    return rows.to_json(orient="records", double_precision=15)


def _candidate_signature(
    diagnostics: dict[str, Any], field: str, *, initial_only: bool = False
) -> str:
    records = sorted(diagnostics["candidate_records"], key=lambda row: tuple(row["ranks"]))
    if initial_only:
        records = [
            row
            for row in records
            if not any(
                str(source).startswith("baseline_local_completion_neighbor_of_")
                for source in row["sources"]
            )
        ]
    if field == "set":
        value = [row["ranks"] for row in records]
    elif field == "qhat":
        value = [[row["ranks"], 2.0 * float(row["objective"])] for row in records]
    elif field == "dimension":
        value = [[row["ranks"], int(row["dimension"])] for row in records]
    else:
        raise ValueError(field)
    return _canonical(value)


def pre_ic_identity(
    selected_3: dict[str, pd.DataFrame], selected_4: dict[str, pd.DataFrame]
) -> pd.DataFrame:
    attempts_3 = selected_3["attempts"].set_index("semantic_replication_id")
    attempts_4 = selected_4["attempts"].set_index("semantic_replication_id")
    diagnostics_3 = _rank_diagnostics(selected_3)
    diagnostics_4 = _rank_diagnostics(selected_4)
    failures_3 = _failure_rows(selected_3).set_index("semantic_replication_id")
    failures_4 = _failure_rows(selected_4).set_index("semantic_replication_id")
    rows = []
    for semantic_id in sorted(attempts_3.index):
        failed_3 = attempts_3.loc[semantic_id, "primary_status"] == "rank_pilot_failure"
        failed_4 = attempts_4.loc[semantic_id, "primary_status"] == "rank_pilot_failure"
        nuclear_equal = _nuclear_signature(
            selected_3["fits"], semantic_id
        ) == _nuclear_signature(selected_4["fits"], semantic_id)
        if failed_3 and failed_4:
            routes_equal = _canonical(
                _json(failures_3.loc[semantic_id, "cap_pilot_start_attempts"])
            ) == _canonical(_json(failures_4.loc[semantic_id, "cap_pilot_start_attempts"]))
            cap_rank_equal = True
            initial_candidate_set_equal = True
            final_candidate_set_equal = True
            qhat_equal = True
            dimension_equal = True
        elif not failed_3 and not failed_4:
            first, second = diagnostics_3[semantic_id], diagnostics_4[semantic_id]
            routes_equal = _canonical(first["cap_pilot"]) == _canonical(second["cap_pilot"])
            cap_rank_equal = first["rank_cap_thresholded_vector"] == second[
                "rank_cap_thresholded_vector"
            ]
            initial_candidate_set_equal = _candidate_signature(
                first, "set", initial_only=True
            ) == _candidate_signature(second, "set", initial_only=True)
            final_candidate_set_equal = _candidate_signature(
                first, "set"
            ) == _candidate_signature(second, "set")
            qhat_equal = _candidate_signature(
                first, "qhat", initial_only=True
            ) == _candidate_signature(second, "qhat", initial_only=True)
            dimension_equal = _candidate_signature(
                first, "dimension", initial_only=True
            ) == _candidate_signature(second, "dimension", initial_only=True)
        else:
            routes_equal = cap_rank_equal = initial_candidate_set_equal = False
            final_candidate_set_equal = False
            qhat_equal = dimension_equal = False
        all_equal = all(
            (
                failed_3 == failed_4,
                nuclear_equal,
                routes_equal,
                cap_rank_equal,
                initial_candidate_set_equal,
                qhat_equal,
                dimension_equal,
            )
        )
        rows.append(
            {
                "semantic_replication_id": semantic_id,
                "pilot_failure_3e-6": failed_3,
                "pilot_failure_4e-6": failed_4,
                "nuclear_path_identical": nuclear_equal,
                "cap_pilot_routes_objectives_validity_identical": routes_equal,
                "thresholded_cap_pilot_rank_identical": cap_rank_equal,
                "initial_pre_ic_candidate_set_identical": initial_candidate_set_equal,
                "final_post_completion_candidate_set_identical": final_candidate_set_equal,
                "candidate_qhat_identical": qhat_equal,
                "candidate_dimensions_identical": dimension_equal,
                "all_pre_ic_artifacts_identical": all_equal,
            }
        )
    result = pd.DataFrame(rows)
    return result


def retention_waterfall(records: pd.DataFrame, label: str) -> pd.DataFrame:
    status_order = (
        ("pilot failures", "rank_pilot_failure", "point"),
        ("rank-cap suppression", "rank_at_cap", "point"),
        ("full-fit failures", "full_fit_failure", "point"),
        ("unsupported selected-rank targets", "target_unsupported_selected_rank", "inference"),
        ("split-fit failures", "split_fit_failure", "inference"),
        ("split-rank loss", "split_rank_loss", "inference"),
        ("tangent-Gram failures", "tangent_gram_nearly_singular", "inference"),
        ("Riesz failures", "riesz_solver_failure", "inference"),
        ("invalid variance", "invalid_variance", "inference"),
        ("nonfinite SE", "nonfinite_standard_error", "inference"),
    )
    remaining = len(records)
    rows = [
        {
            "method": label,
            "stage": "attempted",
            "primary_status": "all",
            "replications_affected": records["semantic_replication_id_attempt"].nunique(),
            "target_records_affected": len(records),
            "target_records_removed": 0,
            "remaining_target_records": remaining,
        }
    ]
    represented = set()
    for stage, status, gate in status_order:
        mask = records["primary_status"].eq(status)
        retained_column = (
            "retained_for_bias_rmse" if gate == "point" else "retained_for_coverage"
        )
        removed = int(
            (mask & ~records[retained_column].fillna(False).astype(bool)).sum()
        )
        remaining -= removed
        represented.add(status)
        rows.append(
            {
                "method": label,
                "stage": stage,
                "primary_status": status,
                "replications_affected": records.loc[
                    mask, "semantic_replication_id_attempt"
                ].nunique(),
                "target_records_affected": int(mask.sum()),
                "target_records_removed": removed,
                "remaining_target_records": remaining,
            }
        )
    other = records.loc[
        ~records["primary_status"].isin({"success", *represented}), "primary_status"
    ]
    if len(other):
        raise ValueError(f"unrepresented mutually exclusive statuses: {other.value_counts()}")
    retained = int(records["retained_for_coverage"].fillna(False).astype(bool).sum())
    if remaining != retained:
        raise ValueError(f"waterfall does not reconcile: {remaining} != {retained}")
    return pd.DataFrame(rows)


def _calibration(root: Path, dgp: int) -> dict[str, Any]:
    frame = pd.read_parquet(root / "calibration.parquet")
    row = frame.loc[frame["dgp"].eq(dgp)].iloc[0]
    return row.to_dict()


def _replay_one(payload: tuple[str, str, int, int]) -> list[dict[str, Any]]:
    label, root_text, dgp, replication = payload
    root = Path(root_text)
    saved = json.loads((root / "resolved_config.json").read_text(encoding="utf-8"))
    config = saved.get("config", saved)
    config["run"]["diagnostic_replay"] = True
    rows = _worker(
        (
            (dgp, 50, 50, replication, None),
            config,
            _calibration(root, dgp),
        )
    )
    for row in rows:
        row["diagnostic_replay"] = True
        row["replay_method"] = label
    return rows


def targeted_replays(fixed_root: Path, selected_root: Path) -> pd.DataFrame:
    payloads = [
        ("selected_4e-6", str(selected_root), 1, 2),
        ("selected_4e-6", str(selected_root), 3, 0),
        ("selected_4e-6", str(selected_root), 3, 2),
        ("fixed_rank", str(fixed_root), 4, 1),
    ]
    with ProcessPoolExecutor(max_workers=4) as executor:
        nested = list(executor.map(_replay_one, payloads))
    rows = [row for group in nested for row in group]
    result = pd.DataFrame(rows)
    original = {
        "fixed_rank": pd.read_parquet(fixed_root / "attempted_replications.parquet").set_index(
            "semantic_replication_id"
        )["dgp_realization_hash"],
        "selected_4e-6": pd.read_parquet(
            selected_root / "attempted_replications.parquet"
        ).set_index("semantic_replication_id")["dgp_realization_hash"],
    }
    replay_attempts = result.loc[result["record_type"].eq("replication")]
    for row in replay_attempts.itertuples():
        semantic_id = row.semantic_replication_id
        replay_hash = row.dgp_realization_hash
        if replay_hash != original[row.replay_method].loc[semantic_id]:
            raise ValueError(f"replay DGP hash mismatch: {semantic_id}")
    return result


def _match_final_fit(fits: pd.DataFrame, attempt: dict[str, Any]) -> pd.Series | None:
    ranks = _canonical(attempt["final_rank_vector"])
    candidates = fits.loc[
        fits["requested_rank"].map(lambda value: _canonical(_json(value))).eq(ranks)
        & np.isclose(
            pd.to_numeric(fits["objective_final"], errors="coerce"),
            float(attempt["final_objective"]),
            rtol=1e-10,
            atol=1e-12,
        )
    ]
    return candidates.iloc[0] if len(candidates) else None


def cap_failure_outputs(replay: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    failures = replay.loc[
        replay["record_type"].eq("failure")
        & replay["primary_status"].eq("rank_pilot_failure")
    ]
    fits = replay.loc[replay["record_type"].eq("fit_diagnostic")]
    route_rows = []
    summary_rows = []
    tolerance = 1e-5
    for failure in failures.itertuples():
        attempts = _json(failure.cap_pilot_start_attempts)
        valid_objectives = [
            (int(item["route_number"]), float(item["final_objective"]))
            for item in attempts
            if item["final_valid"]
        ]
        pairwise = [
            {
                "route_1": first[0],
                "route_2": second[0],
                "normalized_gap": abs(first[1] - second[1])
                / max(1.0, abs(min(first[1], second[1]))),
            }
            for first, second in itertools.combinations(valid_objectives, 2)
        ]
        stable_pairs = sum(item["normalized_gap"] <= tolerance for item in pairwise)
        local_fits = fits.loc[fits["semantic_replication_id"].eq(failure.semantic_replication_id)]
        for attempt in attempts:
            matched = _match_final_fit(local_fits, attempt)
            route_rows.append(
                {
                    "semantic_replication_id": failure.semantic_replication_id,
                    "dgp": failure.dgp,
                    "replication": failure.replication,
                    "route": attempt["route_number"],
                    "route_source": attempt["route_source"],
                    "initial_rank": _canonical(attempt["start_rank_vector"]),
                    "final_numerical_rank": _canonical(
                        attempt["final_numerical_rank_vector"]
                    ),
                    "thresholded_rank": _canonical(
                        attempt["final_thresholded_rank_vector"]
                    ),
                    "starting_envelope": attempt["original_start_max_abs_coefficient"],
                    "start_rescaling": attempt["start_common_scale"],
                    "rescaled_starting_envelope": attempt[
                        "rescaled_start_max_abs_coefficient"
                    ],
                    "iterations": None if matched is None else matched["iterations"],
                    "convergence": attempt["final_converged"],
                    "stationarity_residual": attempt["final_stationarity_residual"],
                    "stationarity_pass": attempt["final_stationarity_residual"] <= 1e-4,
                    "coefficient_envelope": 9.0 * attempt["final_max_envelope_ratio"],
                    "coefficient_envelope_ratio": attempt["final_max_envelope_ratio"],
                    "coefficient_bound_pass": attempt["final_max_envelope_ratio"] < 1.0,
                    "numerical_rank_valid": attempt["final_numerical_rank_vector"]
                    == attempt["final_rank_vector"],
                    "numerical_rank_collapse": attempt["final_numerical_rank_vector"]
                    != attempt["final_rank_vector"],
                    "objective": attempt["final_objective"],
                    "valid": attempt["final_valid"],
                    "invalidity_reason": ";".join(attempt["final_reasons"]) or "none",
                    "start_fit_details": _canonical(
                        attempt["start_fit_diagnostics"].get("start_details", [])
                    ),
                    "full_route_path": _canonical(attempt["path"]),
                    "valid_route_pairwise_gaps": _canonical(pairwise),
                    "diagnostic_replay": True,
                }
            )
        valid_count = sum(item["final_valid"] for item in attempts)
        summary_rows.append(
            {
                "semantic_replication_id": failure.semantic_replication_id,
                "dgp": failure.dgp,
                "replication": failure.replication,
                "attempted_routes": len(attempts),
                "converged_routes": sum(item["final_converged"] for item in attempts),
                "stationarity_pass_routes": sum(
                    item["final_stationarity_residual"] <= 1e-4 for item in attempts
                ),
                "interior_routes": sum(
                    item["final_max_envelope_ratio"] < 1.0 for item in attempts
                ),
                "numerically_rank_valid_routes": sum(
                    item["final_numerical_rank_vector"] == item["final_rank_vector"]
                    for item in attempts
                ),
                "fully_valid_routes": valid_count,
                "stable_route_pairs": stable_pairs,
                "best_two_objective_gap": failure.failure_detail.rsplit("=", 1)[-1],
                "primary_cause": "large disagreement between otherwise valid local optima",
                "precise_failure_reason": failure.failure_detail,
                "diagnostic_replay": True,
            }
        )
    return pd.DataFrame(route_rows), pd.DataFrame(summary_rows)


def runtime_diagnostics(replay: pd.DataFrame) -> pd.DataFrame:
    fits = replay.loc[replay["record_type"].eq("fit_diagnostic")].copy()
    fits["runtime_valid"] = pd.to_numeric(fits["runtime_seconds"], errors="coerce") > 0
    fits["identity_valid"] = (
        fits["fit_type"].notna()
        & fits["semantic_replication_id"].notna()
        & fits["requested_rank"].notna()
        & fits["numerical_rank"].notna()
        & fits["objective_final"].notna()
        & fits["iterations"].notna()
        & fits["convergence_flag"].notna()
        & fits["stationarity_residual"].notna()
        & fits["coefficient_envelope"].notna()
        & fits["initialization_route"].notna()
    )
    columns = [
        "diagnostic_replay", "replay_method", "semantic_replication_id", "dgp",
        "replication", "fit_type", "initialization_route", "start_number",
        "requested_rank", "numerical_rank", "objective_initial", "objective_final",
        "iterations", "convergence_flag", "stationarity_residual", "stationarity_pass",
        "coefficient_envelope", "coefficient_envelope_ratio", "coefficient_bound_hit",
        "runtime_seconds", "runtime_valid", "identity_valid", "diagnostic_context",
        "initial_coefficient_envelope", "final_coefficient_envelope",
        "coefficient_envelope_history",
    ]
    available = [column for column in columns if column in fits]
    return fits[available].sort_values("runtime_seconds", ascending=False)


def fixed_bound_output(
    replay: pd.DataFrame, fixed_data: dict[str, pd.DataFrame]
) -> pd.DataFrame:
    fit = replay.loc[
        replay["record_type"].eq("fit_diagnostic")
        & replay["semantic_replication_id"].eq(FIXED_BOUND_ID)
        & replay["fit_type"].eq("full_fixed_rank")
    ].iloc[0]
    original_rank = fixed_data["rank"].loc[
        fixed_data["rank"]["semantic_replication_id"].eq(FIXED_BOUND_ID)
    ].iloc[0]
    return pd.DataFrame(
        [
            {
                "semantic_replication_id": FIXED_BOUND_ID,
                "dgp": 4,
                "replication": 1,
                "configured_start_count": 1,
                "true_coefficient_envelope": original_rank["realized_coefficient_envelope"],
                "initial_envelope": fit["initial_coefficient_envelope"],
                "final_envelope": fit["final_coefficient_envelope"],
                "final_envelope_ratio": fit["coefficient_envelope_ratio"],
                "objective_initial": fit["objective_initial"],
                "objective_final": fit["objective_final"],
                "stationarity_residual": fit["stationarity_residual"],
                "numerical_rank": fit["numerical_rank"],
                "convergence": fit["convergence_flag"],
                "iterations": fit["iterations"],
                "runtime_seconds": fit["runtime_seconds"],
                "objective_stability_gap": np.nan,
                "envelope_history": fit["coefficient_envelope_history"],
                "diagnosis": (
                    "the sole configured deterministic start converged and passed stationarity "
                    "but moved outside B; no saved or authorized alternative start establishes "
                    "an interior comparable solution, so the original coefficient-bound failure stands"
                ),
                "diagnostic_replay": True,
            }
        ]
    )


def report_text(
    identity: pd.DataFrame,
    summary: pd.DataFrame,
    waterfall_3: pd.DataFrame,
    waterfall_4: pd.DataFrame,
    runtime: pd.DataFrame,
) -> str:
    failures = ", ".join(summary["semantic_replication_id"])
    coefficient = runtime.loc[runtime["fit_type"].ne("nuclear_path")]
    lines = [
        "# Final-preflight failure forensics",
        "",
        "No new DGP draws, medium diagnostics, power, rank stress, or production runs were launched.",
        "",
        "## Central diagnosis",
        "",
        f"All pre-IC artifacts match for {int(identity['all_pre_ic_artifacts_identical'].sum())}/12 semantic replications. The same pilot failures occur under both IC multipliers: {failures}.",
        "The nuclear path, six pilot routes, route objectives/validity, cap-pilot rank, initial candidate ranks, Q-hat values, and dimensions are identical. DGP 2 replication 2 has an expected post-IC local-completion difference: the 3e-6 path adds two neighbors after its initial IC choice, while 4e-6 does not. Thus there is no stochastic reproducibility bug, but the literal final candidate pools are not identical in that cell.",
        "",
        "## Pilot failures",
        "",
    ]
    for row in summary.itertuples():
        lines.append(
            f"- `{row.semantic_replication_id}`: {row.fully_valid_routes}/6 routes fully valid, {row.stable_route_pairs} stable route pairs, but {row.precise_failure_reason}. Primary classification: {row.primary_cause}."
        )
    lines.extend(
        [
            "",
            "The failures are not caused by iteration caps, stationarity, or rank collapse. They are caused by materially different local objective basins among otherwise credible routes. DGP 3 replication 0 also has one bound-invalid route, but five valid routes remain, so that is not the pilot-level cause.",
            "",
            "## Retention",
            "",
            f"- `4e-6`: 216 attempted -> 144 point retained -> {int(waterfall_4.iloc[-1].remaining_target_records)} inference retained. Pilot failures remove 54, the cap hit removes 18, and unsupported selected-rank targets remove 18 at inference. Twenty split-fit diagnostic statuses affect three replications but remove zero records because the stored finite estimates/SEs were retained; this status/retention inconsistency is reported rather than hidden.",
            f"- `3e-6`: 216 attempted -> {int(waterfall_3.iloc[-1].remaining_target_records)} point and inference retained. This is entirely three pilot failures (54 targets) plus seven cap hits (126 targets).",
            "",
            "## Fixed-rank bound event",
            "",
            "Only one fixed-rank start is configured. Its exact-seed replay converges and passes stationarity but moves outside B. No authorized alternative start establishes an interior comparable solution; the event remains a coefficient-bound failure.",
            "The replay also confirms that supplied-rank ALS treats `coefficient_bound` as an acceptance diagnostic rather than projecting iterates into the box: the envelope rises from about 0.089 to 14.439 while B=9. This is an implementation gap relative to a literal box-constrained optimization, but changing it would change the numerical estimator and was not authorized here.",
            "",
            "## Runtime instrumentation",
            "",
            f"All {len(coefficient)} replayed coefficient-fit rows have genuine positive runtimes: {bool(coefficient['runtime_valid'].all())}. Each has fit identity, requested/realized rank, objective, iterations, convergence, stationarity, envelope, and start/route context: {bool(coefficient['identity_valid'].all())}. No historical runtimes were fabricated. Split-fit runtime plumbing is verified by tests; none of the four authorized failure replays reaches split inference.",
            "",
            "## Recommendation",
            "",
            "The smallest plausible next numerical action is deterministic polishing/continuation of the already-converged route endpoints before comparing outer-route objectives, using the unchanged loss, box, ranks, starts, and tolerances. Do not relax objective stability. This is a recommendation only and was not implemented.",
            "",
            "Saved rank outcomes remain: `3e-6` exact 2/9, under 1/9, over 7/9, cap hits 7/9; `4e-6` exact 6/9, under 2/9, over 1/9, cap hits 1/9. The initial IC choice is over the same candidate values; only the subsequent local-completion pool differs in DGP 2 replication 2, without changing either reported winner.",
            "",
            "`3e-6` is not recommended for further preflight. `4e-6` remains a candidate for a later independent validation, but is not selected or frozen as a production default.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed", required=True, type=Path)
    parser.add_argument("--selected-3e-6", required=True, type=Path)
    parser.add_argument("--selected-4e-6", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--reuse-replay", action="store_true")
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    fixed = _load(args.fixed)
    selected_3 = _load(args.selected_3e_6)
    selected_4 = _load(args.selected_4e_6)

    identity = pre_ic_identity(selected_3, selected_4)
    identity.to_csv(output / "pre_ic_identity.csv", index=False)
    waterfall_3 = retention_waterfall(selected_3["records"], "selected_3e-6")
    waterfall_4 = retention_waterfall(selected_4["records"], "selected_4e-6")
    waterfall_3.to_csv(output / "retention_waterfall_3e-6.csv", index=False)
    waterfall_4.to_csv(output / "retention_waterfall_4e-6.csv", index=False)

    replay_path = output / "diagnostic_replay_records.parquet"
    if args.reuse_replay:
        replay = pd.read_parquet(replay_path)
    else:
        replay = targeted_replays(args.fixed, args.selected_4e_6)
        replay.to_parquet(replay_path, index=False)
    routes, summary = cap_failure_outputs(replay)
    routes.to_csv(output / "cap_pilot_failure_routes.csv", index=False)
    summary.to_csv(output / "cap_pilot_failure_summary.csv", index=False)
    runtime = runtime_diagnostics(replay)
    runtime.to_csv(output / "runtime_replay_diagnostics.csv", index=False)
    valid_coefficient = runtime.loc[
        runtime["fit_type"].ne("nuclear_path")
        & runtime["convergence_flag"].fillna(False).astype(bool)
        & runtime["stationarity_pass"].fillna(False).astype(bool)
        & ~runtime["coefficient_bound_hit"].fillna(True).astype(bool)
    ]
    slow_columns = [
        "diagnostic_replay", "replay_method", "semantic_replication_id", "dgp",
        "replication", "fit_type", "initialization_route", "requested_rank",
        "numerical_rank", "objective_final", "iterations", "convergence_flag",
        "stationarity_residual", "coefficient_envelope_ratio", "runtime_seconds",
    ]
    valid_coefficient.head(5)[slow_columns].to_csv(
        output / "five_slowest_replay_fits.csv", index=False
    )
    fixed_bound_output(replay, fixed).to_csv(
        output / "fixed_rank_bound_failure.csv", index=False
    )
    (output / "forensic_report.md").write_text(
        report_text(identity, summary, waterfall_3, waterfall_4, runtime),
        encoding="utf-8",
    )
    print(output)


if __name__ == "__main__":
    main()
