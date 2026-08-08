"""Write nonnumerical method-comparison schemas and CLI documentation."""

from __future__ import annotations

from pathlib import Path

from dynamic_panel_econ.cli import write_cli_help
from dynamic_panel_econ.mc_accounting import write_schema_files


def main() -> None:
    root = Path("results/mc/method_comparison")
    paths = write_schema_files(root)
    help_path = root / "cli_help.txt"
    write_cli_help(help_path)
    readme = root / "README.md"
    readme.write_text(
        "# Method-comparison reporting namespace\n\n"
        "This namespace contains schemas, interface documentation, filtering/reconciliation "
        "audits, and eventual fixed-versus-selected reporting artifacts. It contains no "
        "fabricated Monte Carlo results.\n",
        encoding="utf-8",
    )
    print("\n".join(str(path) for path in [*paths, help_path, readme]))


if __name__ == "__main__":
    main()
