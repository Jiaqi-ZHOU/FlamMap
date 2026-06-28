#!/usr/bin/env python
"""Sweep the CAFT threshold over a range and tabulate LFL/UFL for every .dat.

Re-extracts the flammability limits from existing CAFT grids at each threshold in
[--lo, --hi] step --step, WITHOUT re-running the pipeline (the .dat grid does not
depend on the threshold; only the final extraction does). Output is a long-format CSV
  stem, formula, threshold_K, LFL_percent, UFL_percent
i.e. one row per (molecule, threshold) -- ready to plot LFL/UFL vs threshold and pick
the cutoff(s) that best match experiment.

The work is embarrassingly parallel across .dat files, so on a single (e.g. debug)
node this saturates all allocated cores -- no need for HyperQueue.

Usage:
  # interactive on a debug node (uses all allocated cores by default):
  python threshold_sweep.py --output-dir <DIR> --out sweep.csv
  # explicit ranges / worker count:
  python threshold_sweep.py --output-dir <DIR> --lo 1100 --hi 1800 --step 50 --jobs 128
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")  # headless: no display needed for contouring

SRC_DIR = Path(__file__).resolve().parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from flammability.pd_analysis import compute_flammability_limits

try:
    from tqdm import tqdm
except ImportError:  # tqdm optional
    def tqdm(it, **kw):
        return it


def _find_dats(output_dir: Path) -> list[Path]:
    dat_dir = output_dir / "dat"
    if dat_dir.is_dir() and any(dat_dir.glob("*.dat")):
        return sorted(dat_dir.glob("*.dat"))
    return sorted(output_dir.glob("*.dat"))


def _formula(output_dir: Path, stem: str) -> str:
    jp = output_dir / "json" / f"{stem}.json"
    if jp.is_file():
        try:
            return json.load(open(jp, encoding="utf-8")).get("formula", "") or ""
        except (OSError, json.JSONDecodeError):
            pass
    return ""


# --- worker (one .dat -> rows for every threshold) ---------------------------------
_CTX = {}


def _init(output_dir_str, thresholds):
    _CTX["output_dir"] = Path(output_dir_str)
    _CTX["thresholds"] = thresholds


def _process(dat_str):
    dat = Path(dat_str)
    output_dir, thresholds = _CTX["output_dir"], _CTX["thresholds"]
    stem = dat.stem
    formula = _formula(output_dir, stem)
    rows = []
    for thr in thresholds:
        try:
            lfl, ufl, _ = compute_flammability_limits(dat, thr)
        except Exception:
            lfl = ufl = float("nan")
        rows.append((stem, formula, f"{thr:g}", f"{lfl:.4f}", f"{ufl:.4f}"))
    return rows


def _default_jobs() -> int:
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output-dir", required=True, help="Pipeline output dir with dat/<name>.dat (and json/).")
    p.add_argument("--out", default="threshold_sweep.csv", help="Output CSV path (default: threshold_sweep.csv).")
    p.add_argument("--lo", type=float, default=1100.0, help="Lowest threshold in K (default 1100).")
    p.add_argument("--hi", type=float, default=1800.0, help="Highest threshold in K (default 1800).")
    p.add_argument("--step", type=float, default=50.0, help="Threshold step in K (default 50).")
    p.add_argument("--jobs", type=int, default=None, help="Worker processes (default: all allocated cores).")
    args = p.parse_args()

    output_dir = Path(args.output_dir)
    if not output_dir.is_dir():
        raise SystemExit(f"--output-dir does not exist: {output_dir}")

    n_steps = int(round((args.hi - args.lo) / args.step))
    thresholds = [args.lo + i * args.step for i in range(n_steps + 1)]

    dats = _find_dats(output_dir)
    if not dats:
        raise SystemExit(f"No .dat files found in {output_dir} (or {output_dir}/dat).")
    jobs = max(1, args.jobs or _default_jobs())
    print(f"{len(dats)} molecules x {len(thresholds)} thresholds "
          f"({thresholds[0]:g}-{thresholds[-1]:g} K step {args.step:g}) "
          f"on {jobs} worker(s) -> {args.out}", flush=True)

    out_path = Path(args.out)
    dat_strs = [str(d) for d in dats]
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["stem", "formula", "threshold_K", "LFL_percent", "UFL_percent"])
        if jobs == 1:
            _init(str(output_dir), thresholds)
            for ds in tqdm(dat_strs, unit="mol"):
                w.writerows(_process(ds))
        else:
            with ProcessPoolExecutor(
                max_workers=jobs, initializer=_init, initargs=(str(output_dir), thresholds)
            ) as ex:
                for rows in tqdm(ex.map(_process, dat_strs, chunksize=8), total=len(dat_strs), unit="mol"):
                    w.writerows(rows)
    print(f"Done -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
