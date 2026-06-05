#!/usr/bin/env python
"""Re-extract LFL/UFL from existing CAFT .dat grids at a new threshold.

The CAFT temperature grid (``OUTPUT_DIR/dat/<name>.dat``) does NOT depend on
the threshold — only the final LFL/UFL extraction does. So to see how the
flammability limits change with a different cutoff temperature, there is no
need to re-run the expensive pipeline (TAE / frequencies / YAML / CAFT grid).
This script just reads each existing .dat and re-runs the extraction step at
``--threshold``.

Three modes:

1. Batch (default) — single process over the whole directory, writes one CSV:
       python reextract.py --output-dir OUTPUT_DIR --threshold 1400
   -> OUTPUT_DIR/REEXTRACT_<threshold>.csv

2. Single (--dat) — re-extract ONE .dat and write a per-molecule JSON. Used by
   the HyperQueue template (one task per .dat):
       python reextract.py --output-dir OUTPUT_DIR --threshold 1400 --dat foo.dat
   -> OUTPUT_DIR/reextract_<threshold>/<name>.json

3. Collect (--collect) — fold the per-molecule JSONs from mode 2 into one CSV:
       python reextract.py --output-dir OUTPUT_DIR --threshold 1400 --collect
   -> OUTPUT_DIR/REEXTRACT_<threshold>.csv

Each row is: name, stem, formula, threshold_K, LFL_percent, UFL_percent — where
``formula`` is read back from ``OUTPUT_DIR/json/<name>.json`` when present
(blank otherwise). The original SUMMARY.csv / pipeline JSONs are untouched.
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

_FIELDS = ("name", "stem", "formula", "threshold_K", "LFL_percent", "UFL_percent")


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


def _row_for(output_dir: Path, dat: Path, threshold: float) -> dict:
    name = dat.stem
    lfl, ufl, _segments = compute_flammability_limits(dat, threshold)
    return {
        "name": name,
        "stem": name.split("_")[-1],
        "formula": _formula_for(output_dir, name),
        "threshold_K": f"{threshold:g}",
        "LFL_percent": _fmt(lfl),
        "UFL_percent": _fmt(ufl),
    }


def _per_mol_dir(output_dir: Path, threshold: float) -> Path:
    return output_dir / f"reextract_{threshold:g}"


def run_single(output_dir: Path, threshold: float, dat: Path) -> None:
    """HQ task mode: re-extract one .dat, write a per-molecule JSON."""
    if not dat.is_file():
        raise SystemExit(f"--dat does not exist: {dat}")
    row = _row_for(output_dir, dat, threshold)
    out_dir = _per_mol_dir(output_dir, threshold)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json = out_dir / f"{row['name']}.json"
    with open(out_json, "w", encoding="utf-8") as handle:
        json.dump(row, handle, indent=2)
    print(
        f"{row['name']}: LFL={row['LFL_percent'] or 'nan'}%  "
        f"UFL={row['UFL_percent'] or 'nan'}%  -> {out_json}"
    )


def run_collect(output_dir: Path, threshold: float) -> None:
    """Fold per-molecule JSONs (from run_single) into one CSV."""
    per_mol_dir = _per_mol_dir(output_dir, threshold)
    jsons = sorted(per_mol_dir.glob("*.json"))
    if not jsons:
        raise SystemExit(f"No per-molecule JSONs found under {per_mol_dir}.")
    out_csv = output_dir / f"REEXTRACT_{threshold:g}.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for jp in jsons:
            with open(jp, encoding="utf-8") as jf:
                writer.writerow(json.load(jf))
    print(f"Collected {len(jsons)} molecule(s) -> {out_csv}")


def run_batch(output_dir: Path, threshold: float, out_csv: Path) -> None:
    """Single-process mode: walk the whole directory, write one CSV."""
    dat_files = _find_dat_files(output_dir)
    if not dat_files:
        raise SystemExit(
            f"No .dat files found under {output_dir} (or {output_dir}/dat)."
        )
    print(f"Re-extracting {len(dat_files)} molecule(s) at {threshold:g} K...")
    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for dat in dat_files:
            row = _row_for(output_dir, dat, threshold)
            writer.writerow(row)
            print(
                f"  {row['name']}: LFL={row['LFL_percent'] or 'nan'}%  "
                f"UFL={row['UFL_percent'] or 'nan'}%"
            )
    print(f"Wrote {out_csv}")


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
        "--dat",
        default=None,
        help="Re-extract a single .dat (HQ task mode); writes a per-molecule JSON.",
    )
    parser.add_argument(
        "--collect",
        action="store_true",
        help="Fold per-molecule JSONs from --dat runs into REEXTRACT_<threshold>.csv.",
    )
    parser.add_argument(
        "--output-csv",
        default=None,
        help="Override the batch-mode CSV path "
        "(default: OUTPUT_DIR/REEXTRACT_<threshold>.csv).",
    )
    args = parser.parse_args()

    if args.threshold <= 0:
        parser.error(f"--threshold must be positive, got {args.threshold!r}.")
    if args.dat and args.collect:
        parser.error("--dat and --collect are mutually exclusive.")

    output_dir = Path(args.output_dir)
    if not output_dir.is_dir():
        parser.error(f"--output-dir does not exist: {output_dir}")

    if args.dat:
        run_single(output_dir, args.threshold, Path(args.dat))
    elif args.collect:
        run_collect(output_dir, args.threshold)
    else:
        out_csv = (
            Path(args.output_csv)
            if args.output_csv
            else output_dir / f"REEXTRACT_{args.threshold:g}.csv"
        )
        run_batch(output_dir, args.threshold, out_csv)


if __name__ == "__main__":
    main()
