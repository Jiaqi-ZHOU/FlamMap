#!/usr/bin/env python
"""Re-extract LFL/UFL from existing CAFT .dat grids at a new threshold.

The CAFT temperature grid (``OUTPUT_DIR/dat/<name>.dat``) does NOT depend on
the threshold — only the final LFL/UFL extraction does. So to see how the
flammability limits change with a different cutoff temperature, there is no
need to re-run the expensive pipeline (TAE / frequencies / YAML / CAFT grid).
This script just reads each existing .dat and re-runs the extraction step at
``--threshold``.

Usage:
    python reextract.py --output-dir OUTPUT_DIR --threshold 1400

It writes ``OUTPUT_DIR/REEXTRACT_<threshold>.csv`` with one row per molecule:
    name, stem, formula, threshold_K, LFL_percent, UFL_percent
where ``formula`` is read back from ``OUTPUT_DIR/json/<name>.json`` when that
file exists (blank otherwise). The original SUMMARY.csv / JSONs are untouched.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


SRC_DIR = Path(__file__).resolve().parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from flammability.pd_analysis import compute_flammability_limits


def _find_dat_files(output_dir: Path) -> list[Path]:
    """Prefer the batch layout (OUTPUT_DIR/dat/), fall back to flat OUTPUT_DIR."""
    dat_dir = output_dir / "dat"
    if dat_dir.is_dir():
        files = sorted(dat_dir.glob("*.dat"))
        if files:
            return files
    return sorted(output_dir.glob("*.dat"))


def _formula_for(output_dir: Path, name: str) -> str:
    """Read the chemical formula back from the per-molecule JSON, if present."""
    json_path = output_dir / "json" / f"{name}.json"
    if not json_path.is_file():
        return ""
    try:
        with open(json_path, encoding="utf-8") as handle:
            return json.load(handle).get("formula", "") or ""
    except (OSError, json.JSONDecodeError):
        return ""


def _fmt(value: float) -> str:
    return "" if value is None or math.isnan(value) else f"{value:.6f}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Pipeline output directory containing dat/<name>.dat (and json/).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        required=True,
        help="CAFT cutoff temperature in K for LFL/UFL extraction.",
    )
    parser.add_argument(
        "--output-csv",
        default=None,
        help="Override the output CSV path (default: "
        "OUTPUT_DIR/REEXTRACT_<threshold>.csv).",
    )
    args = parser.parse_args()

    if args.threshold <= 0:
        parser.error(f"--threshold must be positive, got {args.threshold!r}.")

    output_dir = Path(args.output_dir)
    if not output_dir.is_dir():
        parser.error(f"--output-dir does not exist: {output_dir}")

    dat_files = _find_dat_files(output_dir)
    if not dat_files:
        parser.error(f"No .dat files found under {output_dir} (or {output_dir}/dat).")

    out_csv = (
        Path(args.output_csv)
        if args.output_csv
        else output_dir / f"REEXTRACT_{args.threshold:g}.csv"
    )

    fields = ("name", "stem", "formula", "threshold_K", "LFL_percent", "UFL_percent")
    print(f"Re-extracting {len(dat_files)} molecule(s) at {args.threshold:g} K...")
    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for dat in dat_files:
            name = dat.stem
            stem = name.split("_")[-1]
            lfl, ufl, _segments = compute_flammability_limits(dat, args.threshold)
            writer.writerow(
                {
                    "name": name,
                    "stem": stem,
                    "formula": _formula_for(output_dir, name),
                    "threshold_K": f"{args.threshold:g}",
                    "LFL_percent": _fmt(lfl),
                    "UFL_percent": _fmt(ufl),
                }
            )
            print(f"  {name}: LFL={_fmt(lfl) or 'nan'}%  UFL={_fmt(ufl) or 'nan'}%")

    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
