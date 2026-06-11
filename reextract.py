#!/usr/bin/env python
"""Re-extract LFL/UFL from existing CAFT .dat grids at a new threshold.

The CAFT temperature grid (``OUTPUT_DIR/dat/<name>.dat``) does NOT depend on
the threshold — only the final LFL/UFL extraction does. So to see how the
flammability limits change with a different cutoff temperature, there is no
need to re-run the expensive pipeline (TAE / frequencies / YAML / CAFT grid).
This script just reads each existing .dat and re-runs the extraction step at
``--threshold`` (or independent ``--lfl-threshold`` / ``--ufl-threshold``).

The output CSV uses the SAME columns as the pipeline's SUMMARY.csv:
    stem, formula, tae_Ha, Hf_298K_kJ, lfl_threshold_K, ufl_threshold_K,
    LFL_percent, UFL_percent, n_freqs, elapsed_s
``LFL_percent`` / ``UFL_percent`` / ``lfl_threshold_K`` / ``ufl_threshold_K``
reflect the new threshold(s); the other columns (tae_Ha, Hf_298K_kJ, n_freqs,
elapsed_s, formula) are read back from ``OUTPUT_DIR/json/<name>.json`` and left
blank if that JSON is absent. The original SUMMARY.csv and pipeline JSONs are
untouched.

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
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


SRC_DIR = Path(__file__).resolve().parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from flammability.pd_analysis import compute_split_flammability_limits

# Same columns and order as batch.py's _SUMMARY_FIELDS, so REEXTRACT_*.csv is a
# drop-in for SUMMARY.csv in downstream tooling.
_FIELDS = (
    "stem",
    "formula",
    "tae_Ha",
    "Hf_298K_kJ",
    "lfl_threshold_K",
    "ufl_threshold_K",
    "LFL_percent",
    "UFL_percent",
    "n_freqs",
    "elapsed_s",
)


def _thr_tag(lfl_threshold: float, ufl_threshold: float) -> str:
    """File/dir tag for a threshold pair: ``1600`` if equal, else ``1600-1250``."""
    if lfl_threshold == ufl_threshold:
        return f"{lfl_threshold:g}"
    return f"{lfl_threshold:g}-{ufl_threshold:g}"


def _find_dat_files(output_dir: Path) -> list[Path]:
    """Prefer the batch layout (OUTPUT_DIR/dat/), fall back to flat OUTPUT_DIR."""
    dat_dir = output_dir / "dat"
    if dat_dir.is_dir():
        files = sorted(dat_dir.glob("*.dat"))
        if files:
            return files
    return sorted(output_dir.glob("*.dat"))


def _pipeline_meta(output_dir: Path, name: str) -> dict:
    """Read the per-molecule pipeline JSON (tae, Hf, freqs, ...) if present."""
    json_path = output_dir / "json" / f"{name}.json"
    if not json_path.is_file():
        return {}
    try:
        with open(json_path, encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}


def _opt_float(value) -> float | str:
    """float(value), or '' for a missing field — keeps blank cells like SUMMARY."""
    return "" if value is None else float(value)


def _row_for(output_dir: Path, dat: Path, lfl_threshold: float, ufl_threshold: float) -> dict:
    """Build a SUMMARY-shaped row: LFL/UFL from the .dat at the new threshold(s),
    everything else carried over from the molecule's pipeline JSON."""
    name = dat.stem
    meta = _pipeline_meta(output_dir, name)
    lfl, ufl, _lfl_segs, _ufl_segs = compute_split_flammability_limits(
        dat, lfl_threshold, ufl_threshold
    )
    freqs = meta.get("freqs")
    return {
        "stem": name.split("_")[-1],
        "formula": meta.get("formula", "") or "",
        "tae_Ha": _opt_float(meta.get("tae")),
        "Hf_298K_kJ": _opt_float(meta.get("Hf_298K_kJ")),
        "lfl_threshold_K": float(lfl_threshold),
        "ufl_threshold_K": float(ufl_threshold),
        "LFL_percent": float(lfl),
        "UFL_percent": float(ufl),
        "n_freqs": len(freqs) if isinstance(freqs, list) else "",
        "elapsed_s": _opt_float(meta.get("elapsed_s")),
    }


def _per_mol_dir(output_dir: Path, lfl_threshold: float, ufl_threshold: float) -> Path:
    return output_dir / f"reextract_{_thr_tag(lfl_threshold, ufl_threshold)}"


def run_single(output_dir: Path, lfl_threshold: float, ufl_threshold: float, dat: Path) -> None:
    """HQ task mode: re-extract one .dat, write a per-molecule JSON."""
    if not dat.is_file():
        raise SystemExit(f"--dat does not exist: {dat}")
    row = _row_for(output_dir, dat, lfl_threshold, ufl_threshold)
    out_dir = _per_mol_dir(output_dir, lfl_threshold, ufl_threshold)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json = out_dir / f"{dat.stem}.json"
    with open(out_json, "w", encoding="utf-8") as handle:
        json.dump(row, handle, indent=2)
    print(
        f"{dat.stem}: LFL={row['LFL_percent']}  UFL={row['UFL_percent']}  "
        f"-> {out_json}"
    )


def run_collect(output_dir: Path, lfl_threshold: float, ufl_threshold: float) -> None:
    """Fold per-molecule JSONs (from run_single) into one SUMMARY-shaped CSV."""
    per_mol_dir = _per_mol_dir(output_dir, lfl_threshold, ufl_threshold)
    jsons = sorted(per_mol_dir.glob("*.json"))
    if not jsons:
        raise SystemExit(f"No per-molecule JSONs found under {per_mol_dir}.")
    out_csv = output_dir / f"REEXTRACT_{_thr_tag(lfl_threshold, ufl_threshold)}.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for jp in jsons:
            with open(jp, encoding="utf-8") as jf:
                writer.writerow(json.load(jf))
    print(f"Collected {len(jsons)} molecule(s) -> {out_csv}")


def run_batch(output_dir: Path, lfl_threshold: float, ufl_threshold: float, out_csv: Path) -> None:
    """Single-process mode: walk the whole directory, write one CSV."""
    dat_files = _find_dat_files(output_dir)
    if not dat_files:
        raise SystemExit(
            f"No .dat files found under {output_dir} (or {output_dir}/dat)."
        )
    if lfl_threshold == ufl_threshold:
        print(f"Re-extracting {len(dat_files)} molecule(s) at {lfl_threshold:g} K...")
    else:
        print(
            f"Re-extracting {len(dat_files)} molecule(s) "
            f"(LFL at {lfl_threshold:g} K, UFL at {ufl_threshold:g} K)..."
        )
    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for dat in dat_files:
            row = _row_for(output_dir, dat, lfl_threshold, ufl_threshold)
            writer.writerow(row)
            print(
                f"  {dat.stem}: LFL={row['LFL_percent']}  UFL={row['UFL_percent']}"
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
        default=None,
        help="CAFT cutoff temperature in K for LFL/UFL extraction. Sets both limits "
        "unless a per-limit override is given. Required unless both --lfl-threshold "
        "and --ufl-threshold are provided.",
    )
    parser.add_argument(
        "--lfl-threshold",
        type=float,
        default=None,
        help="CAFT cutoff in K for the LOWER limit only (default: --threshold).",
    )
    parser.add_argument(
        "--ufl-threshold",
        type=float,
        default=None,
        help="CAFT cutoff in K for the UPPER limit only (default: --threshold).",
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

    # Resolve per-limit thresholds: each falls back to the shared --threshold.
    lfl_threshold = args.lfl_threshold if args.lfl_threshold is not None else args.threshold
    ufl_threshold = args.ufl_threshold if args.ufl_threshold is not None else args.threshold
    if lfl_threshold is None or ufl_threshold is None:
        parser.error(
            "provide --threshold, or both --lfl-threshold and --ufl-threshold."
        )
    for label, value in (("--lfl-threshold", lfl_threshold), ("--ufl-threshold", ufl_threshold)):
        if value <= 0:
            parser.error(f"{label} must be positive, got {value!r}.")
    if args.dat and args.collect:
        parser.error("--dat and --collect are mutually exclusive.")

    output_dir = Path(args.output_dir)
    if not output_dir.is_dir():
        parser.error(f"--output-dir does not exist: {output_dir}")

    if args.dat:
        run_single(output_dir, lfl_threshold, ufl_threshold, Path(args.dat))
    elif args.collect:
        run_collect(output_dir, lfl_threshold, ufl_threshold)
    else:
        out_csv = (
            Path(args.output_csv)
            if args.output_csv
            else output_dir / f"REEXTRACT_{_thr_tag(lfl_threshold, ufl_threshold)}.csv"
        )
        run_batch(output_dir, lfl_threshold, ufl_threshold, out_csv)


if __name__ == "__main__":
    main()
