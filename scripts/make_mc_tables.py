"""Generate deterministic paper-facing table fragments from aggregate files."""

from __future__ import annotations

import argparse
from pathlib import Path

from dynamic_panel_econ.config import config_hash, load_config
from dynamic_panel_econ.monte_carlo import resolve_group_gap
from dynamic_panel_econ.reporting import make_tables


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-root")
    args = parser.parse_args()
    config, _ = resolve_group_gap(load_config(args.config))
    root = Path(args.run_root) if args.run_root else Path(config["run"]["output_root"]) / config["run"]["name"] / config_hash(config)
    for path in make_tables(root):
        print(path)


if __name__ == "__main__":
    main()
