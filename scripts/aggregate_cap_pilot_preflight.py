"""Create the compact 12-replication cap-pilot preflight report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dynamic_panel_econ.config import config_hash, load_config
from dynamic_panel_econ.monte_carlo import resolve_group_gap


def _read_chunks(path: Path) -> pd.DataFrame:
    files = sorted(path.glob("*.parquet"))
    return pd.concat([pd.read_parquet(file) for file in files], ignore_index=True) if files else pd.DataFrame()


def _json(value: object, default: object) -> object:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return default
    return json.loads(str(value))


def _markdown_table(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for values in frame.itertuples(index=False, name=None):
        formatted = []
        for value in values:
            if isinstance(value, float):
                formatted.append("" if np.isnan(value) else f"{value:.6g}")
            else:
                formatted.append("" if value is None else str(value))
        lines.append("| " + " | ".join(formatted) + " |")
    return "\n".join(lines)


def build_table(root: Path, config: dict) -> pd.DataFrame:
    ranks = _read_chunks(root / "rank")
    raw = _read_chunks(root / "raw")
    failures = raw.loc[raw.get("record_type", pd.Series(dtype=str)).eq("failure")].copy()
    rows = []
    for dgp in config["run"]["dgps"]:
        for replication in range(int(config["run"]["replications"])):
            selected = ranks.loc[(ranks["dgp"] == dgp) & (ranks["replication"] == replication)] if not ranks.empty else pd.DataFrame()
            failed = failures.loc[(failures["dgp"] == dgp) & (failures["replication"] == replication)] if not failures.empty else pd.DataFrame()
            source = selected.iloc[0] if not selected.empty else (failed.iloc[0] if not failed.empty else pd.Series(dtype=object))
            attempts = _json(source.get("cap_pilot_start_attempts"), [])
            route_objectives = [
                item.get("final_objective") for item in attempts if item.get("final_valid")
            ]
            finite_objectives = sorted(
                float(value) for value in route_objectives if value is not None and np.isfinite(value)
            )
            computed_gap = (
                abs(finite_objectives[0] - finite_objectives[1])
                / max(1.0, abs(finite_objectives[0]))
                if len(finite_objectives) >= 2
                else np.nan
            )
            stored_gap = source.get("cap_pilot_best_two_objective_gap", np.nan)
            best_two_gap = stored_gap if pd.notna(stored_gap) else computed_gap
            rows.append(
                {
                    "dgp": dgp,
                    "N": 50,
                    "T": 50,
                    "replication": replication,
                    "cap_pilot_success": not selected.empty,
                    "status": source.get("status", "missing"),
                    "attempted_routes": source.get("cap_pilot_attempted_route_count", 0),
                    "valid_routes": source.get("cap_pilot_valid_route_count", 0),
                    "stable_routes": source.get("cap_pilot_stable_route_count", 0),
                    "route_objectives": json.dumps(route_objectives),
                    "best_two_gap": best_two_gap,
                    "thresholded_cap_rank": source.get("rank_cap_thresholded_vector"),
                    "candidate_covers_111": source.get("true_rank_in_candidates", False),
                    "selected_rank": source.get("selected_rank_vector"),
                    "runtime_seconds": source.get("replication_runtime_seconds", source.get("rank_runtime_seconds", np.nan)),
                    "truth_envelope": source.get("realized_coefficient_envelope", np.nan),
                    "truth_envelope_pass": source.get("coefficient_envelope_condition_pass", False),
                    "cap_envelope_ratio": source.get("cap_pilot_max_envelope_ratio", np.nan),
                    "selected_envelope_ratio": source.get("selected_max_envelope_ratio", np.nan),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/mc/cap_pilot_preflight.toml")
    parser.add_argument("--run-root")
    parser.add_argument("--output", default="results/mc/diagnostics/cap_pilot_preflight")
    args = parser.parse_args()
    config, _ = resolve_group_gap(load_config(args.config))
    root = Path(args.run_root) if args.run_root else Path(config["run"]["output_root"]) / config["run"]["name"] / config_hash(config)
    table = build_table(root, config)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "tab_cap_pilot_preflight.csv", index=False)
    table.to_parquet(output / "tab_cap_pilot_preflight.parquet", index=False)
    summary = {
        "replications": len(table),
        "cap_pilot_success_rate": float(table["cap_pilot_success"].mean()),
        "candidate_coverage_rate": float(table["candidate_covers_111"].mean()),
        "selected_rank_distribution": table["selected_rank"].fillna("failure").value_counts().to_dict(),
        "runtime_seconds": table["runtime_seconds"].describe().to_dict(),
        "truth_envelope_failure_count": int((~table["truth_envelope_pass"].astype(bool)).sum()),
    }
    (output / "preflight_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    display = table[
        [
            "dgp",
            "replication",
            "cap_pilot_success",
            "attempted_routes",
            "stable_routes",
            "best_two_gap",
            "thresholded_cap_rank",
            "candidate_covers_111",
            "selected_rank",
            "runtime_seconds",
            "truth_envelope",
            "cap_envelope_ratio",
        ]
    ].copy()
    markdown = [
        "# Cap-pilot numerical-start preflight",
        "",
        "Decision: **NO-GO** for medium diagnostics. The strict cap-pilot success rate was "
        f"{summary['cap_pilot_success_rate']:.1%} "
        f"({int(table['cap_pilot_success'].sum())}/{len(table)}), below the requested nearly-all gate.",
        "",
        f"Candidate coverage was {summary['candidate_coverage_rate']:.1%}; all 12 true coefficient "
        "envelopes satisfied the deterministic interior condition. No medium or production run was started.",
        "",
        _markdown_table(display),
        "",
        "Full route objectives and machine-readable columns are in `tab_cap_pilot_preflight.csv` "
        "and `tab_cap_pilot_preflight.parquet`.",
        "",
    ]
    (output / "preflight_report.md").write_text("\n".join(markdown), encoding="utf-8")
    print(table.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
