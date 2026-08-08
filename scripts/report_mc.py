"""Create reconciled method-comparison tables and figures."""

from __future__ import annotations

import argparse

from dynamic_panel_econ.method_reporting import (
    FIGURE_NAMES,
    TABLE_NAMES,
    report_method_comparison,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-run", action="append", required=True)
    parser.add_argument("--output-dir", default="results/mc/method_comparison")
    parser.add_argument("--tables", nargs="+", default=[], choices=(*TABLE_NAMES, "all"))
    parser.add_argument("--figures", nargs="+", default=[], choices=(*FIGURE_NAMES, "all"))
    parser.add_argument("--table", action="append", choices=(*TABLE_NAMES, "failures"))
    parser.add_argument("--figure", action="append", choices=(*FIGURE_NAMES, "power"))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    tables = args.table or args.tables
    figures = args.figure or args.figures
    tables = ["failure_accounting" if name == "failures" else name for name in tables]
    figures = [name for item in figures for name in (("power_A", "power_B") if item == "power" else (item,))]
    if not tables and not figures:
        tables, figures = ["all"], ["all"]
    result = report_method_comparison(
        args.input_run,
        args.output_dir,
        tables=tables,
        figures=figures,
        overwrite=args.overwrite,
    )
    print(result)


if __name__ == "__main__":
    main()
