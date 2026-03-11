from __future__ import annotations

import argparse

from .config import load_config, validate_config
from .stages import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the single-species flammability pipeline.")
    parser.add_argument("config", help="Path to the pipeline YAML config.")
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate the config and exit without running the pipeline.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    cfg = load_config(args.config)

    errors = validate_config(cfg)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        raise SystemExit(1)

    if args.validate_only:
        print("Config is valid.")
        return

    run_pipeline(cfg)


if __name__ == "__main__":
    main()
