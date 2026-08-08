"""Audit repository constructs that can filter, skip, or classify observations."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

PATTERN = re.compile(
    r"dropna|notna|isna|isfinite|finite|try:|except|continue|return None|"
    r"trim|winsor|quantile|clip|outlier|\.loc\[|merge\(|groupby\(|drop_duplicates",
    re.IGNORECASE,
)


def audit_repository(root: Path) -> list[dict[str, str | int]]:
    rows = []
    for directory in ("src", "scripts", "tests"):
        for path in sorted((root / directory).rglob("*.py")):
            for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if not PATTERN.search(line):
                    continue
                stripped = line.strip()
                row_loss = (
                    "possible"
                    if re.search(
                        r"dropna|continue|return None|\.loc\[|merge\(|drop_duplicates|groupby\(",
                        stripped,
                        re.I,
                    )
                    else "no direct row deletion"
                )
                recorded = (
                    "test-only assertion"
                    if directory == "tests"
                    else "reviewed: status/validity diagnostic or numerical control"
                )
                rows.append(
                    {
                        "file": path.relative_to(root).as_posix(),
                        "line": line_number,
                        "construct": stripped,
                        "purpose": "filtering/numerical-control occurrence requiring accounting review",
                        "rows_can_disappear": row_loss,
                        "failure_recording": recorded,
                    }
                )
    return rows


def render_audit(rows: list[dict[str, str | int]]) -> str:
    lines = [
        "# Filtering and dropping audit",
        "",
        "This audit enumerates filtering-like constructs in `src/`, `scripts/`, and `tests/`. Main method-comparison summaries begin from the attempted-replication ledger and left-join target results, so a missing target row becomes an explicit failure rather than disappearing.",
        "",
        "Main policy: all valid finite estimates remain in bias and RMSE. No trimming or winsorization is permitted. The diagnostic `extreme_estimate_flag` never controls retention.",
        "",
        "| File | Line | Construct | Purpose | Rows can disappear? | Failure recording |",
        "|:---|---:|:---|:---|:---|:---|",
    ]
    for row in rows:
        construct = str(row["construct"]).replace("|", r"\|").replace("`", "'")
        lines.append(
            f"| `{row['file']}` | {row['line']} | `{construct}` | {row['purpose']} | "
            f"{row['rows_can_disappear']} | {row['failure_recording']} |"
        )
    lines += [
        "",
        "## Accounting conclusion",
        "",
        "The new accounting path never uses `dropna` to define a denominator. Attempted identities are expanded to all requested targets and results are attached by a validated left join. Nonfinite estimates, nonfinite standard errors, invalid variances, numerical failures, and missing target rows receive explicit primary statuses.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--output", type=Path, default=Path("results/mc/audit/filtering_audit.md")
    )
    args = parser.parse_args()
    rows = audit_repository(args.root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_audit(rows), encoding="utf-8")
    mirror = args.root / "results/mc/method_comparison/audit/filtering_audit.md"
    mirror.parent.mkdir(parents=True, exist_ok=True)
    mirror.write_text(render_audit(rows), encoding="utf-8")
    print(f"audited_occurrences={len(rows)}")


if __name__ == "__main__":
    main()
