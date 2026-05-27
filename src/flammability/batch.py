"""Batch driver: run the flammability pipeline over many XYZ files.

Per-molecule outputs (yaml / dat / pdf / json) keep their current shape and
land in ``output_dir`` named by the xyz filename stem (not the chemical
formula, which can collide across isomers in a batch). The batch driver
additionally writes:

- ``output_dir/_summary.csv`` — one row per successful molecule.
- ``output_dir/_failed.csv`` — one row per failed molecule (stem, error_type,
  error, elapsed_s). Use single-molecule mode (`run.py <one.xyz>`) to get a
  full traceback when debugging a specific failure.

The ``stem`` column is the trailing ``_``-separated segment of the xyz
filename (e.g. ``"C2H2_ca2cc2cc"`` → ``"ca2cc2cc"``) — typically the unique
hash that distinguishes isomers in a dataset.

A failure is logged and the loop continues; the whole batch never aborts
because of a single bad input. Re-running with ``--skip-existing`` resumes
from wherever the previous attempt died.
"""

from __future__ import annotations

import csv
import json
import multiprocessing as mp
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .config import DEFAULT_OUTPUT_DIR, build_config, validate_config


_SUMMARY_FIELDS = (
    "stem",
    "formula",
    "tae_Ha",
    "Hf_298K_kJ",
    "LFL_percent",
    "UFL_percent",
    "n_freqs",
    "elapsed_s",
)

_FAILED_FIELDS = (
    "stem",
    "error_type",
    "error",
    "elapsed_s",
)


def _read_input_list(path: Path) -> list[Path]:
    """Read XYZ paths, one per line. Relative paths are resolved against the
    input list file's parent directory (not against the process cwd), so the
    list is portable regardless of where the batch is launched from."""
    base = path.parent
    xyz_paths: list[Path] = []
    with path.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            p = Path(line).expanduser()
            if not p.is_absolute():
                p = (base / p).resolve()
            xyz_paths.append(p)
    return xyz_paths


def _worker_init(threads_per_worker: int) -> None:
    # Cap threading inside each worker so N workers × M threads don't oversubscribe.
    # MUST be set before numpy / torch / pyscf are first imported.
    os.environ["OMP_NUM_THREADS"] = str(threads_per_worker)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads_per_worker)
    os.environ["MKL_NUM_THREADS"] = str(threads_per_worker)
    os.environ["NUMEXPR_NUM_THREADS"] = str(threads_per_worker)


def _run_one(
    *,
    xyz: str,
    mode: str,
    case_yaml: str | None,
    output_dir: str,
    plot_map: bool,
    hip_checkpoint: str | None,
    hip_device: str,
    skala_device: str,
) -> dict:
    """Run the pipeline for one molecule; never raise."""
    from .stages import run_pipeline

    xyz_path = Path(xyz)
    name = xyz_path.stem                  # full file basename, used for output files + console
    stem = name.split("_")[-1]            # trailing hash segment, used as CSV identifier
    start = time.perf_counter()

    try:
        cfg = build_config(
            xyz_geom=xyz_path,
            mode=mode,
            case_yaml=case_yaml,
            output_dir=output_dir,
            plot_map=plot_map,
            hip_checkpoint=hip_checkpoint,
            hip_device=hip_device,
            skala_device=skala_device,
        )
        errors = validate_config(cfg)
        if errors:
            raise ValueError("; ".join(errors))

        summary = run_pipeline(cfg, quiet=True, case_name_override=name)
        elapsed = time.perf_counter() - start

        return {
            "status": "ok",
            "name": name,    # for console display; dropped from CSV via extrasaction="ignore"
            "stem": stem,
            "formula": summary["formula"],
            "tae_Ha": float(summary["tae"]),
            "Hf_298K_kJ": float(summary["Hf_298K_kJ"]),
            "LFL_percent": float(summary["LFL_percent"]),
            "UFL_percent": float(summary["UFL_percent"]),
            "n_freqs": len(summary.get("freqs") or []),
            "elapsed_s": round(elapsed, 3),
        }
    except Exception as exc:
        # Keep one-line error suitable for CSV; full traceback is reproducible
        # by re-running the offending molecule with `run.py <xyz>` standalone.
        return {
            "status": "fail",
            "name": name,
            "stem": stem,
            "error_type": type(exc).__name__,
            "error": " ".join(str(exc).split()),  # collapse internal newlines/whitespace
            "elapsed_s": round(time.perf_counter() - start, 3),
        }


def run_batch(
    *,
    input_list: str,
    mode: str,
    case_yaml: str | None,
    output_dir: str | None,
    plot_map: bool,
    hip_checkpoint: str | None,
    hip_device: str,
    skala_device: str,
    jobs: int,
    skip_existing: bool,
    validate_only: bool,
) -> None:
    input_list_path = Path(input_list).expanduser().resolve()
    if not input_list_path.is_file():
        print(f"ERROR: --input-list not found: {input_list_path}")
        raise SystemExit(1)

    out_dir = (
        Path(output_dir).expanduser().resolve()
        if output_dir is not None
        else DEFAULT_OUTPUT_DIR
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    xyz_paths = _read_input_list(input_list_path)
    if not xyz_paths:
        print(f"ERROR: --input-list is empty: {input_list_path}")
        raise SystemExit(1)

    if skip_existing:
        before = len(xyz_paths)
        xyz_paths = [p for p in xyz_paths if not (out_dir / "json" / f"{p.stem}.json").is_file()]
        skipped = before - len(xyz_paths)
        if skipped:
            print(f"--skip-existing: skipping {skipped} input(s) that already have output JSON.")

    if validate_only:
        bad = [p for p in xyz_paths if not p.is_file()]
        if bad:
            for p in bad:
                print(f"ERROR: missing XYZ: {p}")
            raise SystemExit(1)
        print(f"Inputs are valid. {len(xyz_paths)} XYZ files queued.")
        return

    n_total = len(xyz_paths)
    # Prefer the CPU set the OS scheduler assigned to this process — under
    # SLURM that's exactly what `--cpus-per-task` allocated, which is what we
    # want to share across workers. `os.cpu_count()` would return the whole
    # node's logical CPUs, double-counting hyperthreads and ignoring cgroup
    # limits, leading to oversubscription on shared nodes.
    try:
        nproc = len(os.sched_getaffinity(0))
    except AttributeError:  # non-Linux fallback
        nproc = os.cpu_count() or 1
    threads_per_worker = max(1, nproc // max(1, jobs))
    print(
        f"Batch: {n_total} molecules, {jobs} worker(s), "
        f"{threads_per_worker} thread(s)/worker (nproc={nproc}), output={out_dir}"
    )

    summary_csv = out_dir / "_summary.csv"
    failed_csv = out_dir / "_failed.csv"
    new_summary = not summary_csv.is_file()
    new_failed = not failed_csv.is_file()

    one_kwargs = dict(
        mode=mode,
        case_yaml=case_yaml,
        output_dir=str(out_dir),
        plot_map=plot_map,
        hip_checkpoint=hip_checkpoint,
        hip_device=hip_device,
        skala_device=skala_device,
    )

    t0 = time.perf_counter()
    n_ok = 0
    n_fail = 0

    with summary_csv.open("a", newline="", encoding="utf-8") as sf, failed_csv.open(
        "a", newline="", encoding="utf-8"
    ) as ff:
        sum_writer = csv.DictWriter(sf, fieldnames=_SUMMARY_FIELDS, extrasaction="ignore")
        fail_writer = csv.DictWriter(ff, fieldnames=_FAILED_FIELDS, extrasaction="ignore")
        if new_summary:
            sum_writer.writeheader()
            sf.flush()
        if new_failed:
            fail_writer.writeheader()
            ff.flush()

        if jobs <= 1:
            # Single-process loop — amortises startup across all molecules in
            # one Python interpreter. Best when jobs=1 (saves spawn overhead).
            for i, xyz in enumerate(xyz_paths, 1):
                result = _run_one(xyz=str(xyz), **one_kwargs)
                _record(result, sum_writer, fail_writer, i, n_total)
                if result["status"] == "ok":
                    n_ok += 1
                else:
                    n_fail += 1
                sf.flush()
                ff.flush()
        else:
            # Spawn so CUDA contexts (if any) don't get inherited across forks.
            ctx = mp.get_context("spawn")
            with ProcessPoolExecutor(
                max_workers=jobs,
                mp_context=ctx,
                initializer=_worker_init,
                initargs=(threads_per_worker,),
            ) as pool:
                futures = {
                    pool.submit(_run_one, xyz=str(p), **one_kwargs): p
                    for p in xyz_paths
                }
                for i, fut in enumerate(as_completed(futures), 1):
                    result = fut.result()
                    _record(result, sum_writer, fail_writer, i, n_total)
                    if result["status"] == "ok":
                        n_ok += 1
                    else:
                        n_fail += 1
                    sf.flush()
                    ff.flush()

    elapsed = time.perf_counter() - t0
    print(
        f"\nBatch complete in {elapsed:.1f} s: {n_ok} ok, {n_fail} failed. "
        f"Summary: {summary_csv}"
    )
    if n_fail:
        print(f"Failures logged to: {failed_csv}")


_CSV_FLOAT_FIELDS_3DP = ("tae_Ha", "Hf_298K_kJ", "LFL_percent", "UFL_percent", "elapsed_s")


def _record(
    result: dict,
    sum_writer: csv.DictWriter,
    fail_writer: csv.DictWriter,
    i: int,
    total: int,
) -> None:
    if result["status"] == "ok":
        row = {
            k: (f"{v:.3f}" if k in _CSV_FLOAT_FIELDS_3DP and isinstance(v, float) else v)
            for k, v in result.items()
        }
        sum_writer.writerow(row)
        print(
            f"[{i}/{total}] {result['name']:24s} OK  "
            f"Hf={result['Hf_298K_kJ']:>8.2f} kJ/mol  "
            f"LFL={result['LFL_percent']:.3f}%  UFL={result['UFL_percent']:.3f}%  "
            f"({result['elapsed_s']:.1f}s)"
        )
    else:
        row = {
            k: (f"{v:.3f}" if k == "elapsed_s" and isinstance(v, float) else v)
            for k, v in result.items()
        }
        fail_writer.writerow(row)
        print(
            f"[{i}/{total}] {result['name']:24s} FAIL  "
            f"{result['error_type']}: {result['error']}  ({result['elapsed_s']:.1f}s)"
        )


def collect_summary(*, input_list: str, output_dir: str | None) -> None:
    """Rebuild ``_summary.csv`` and ``_failed.csv`` from per-molecule JSONs.

    Use after a HyperQueue (or other external) batch run that produced
    ``<output_dir>/json/<name>.json`` files but no aggregated CSV. Molecules
    listed in ``input_list`` but missing a JSON are treated as failures.
    The ``elapsed_s`` column is left blank because external runners don't
    feed back per-task timing — check HQ's task report for that.
    """
    out_dir = (
        Path(output_dir).expanduser().resolve()
        if output_dir is not None
        else DEFAULT_OUTPUT_DIR
    )
    json_dir = out_dir / "json"
    if not json_dir.is_dir():
        print(f"ERROR: no json/ subdir found in {out_dir}")
        raise SystemExit(1)

    input_list_path = Path(input_list).expanduser().resolve()
    if not input_list_path.is_file():
        print(f"ERROR: --collect-summary list not found: {input_list_path}")
        raise SystemExit(1)

    xyz_paths = _read_input_list(input_list_path)
    summary_csv = out_dir / "_summary.csv"
    failed_csv = out_dir / "_failed.csv"

    n_ok = 0
    n_fail = 0
    with summary_csv.open("w", newline="", encoding="utf-8") as sf, failed_csv.open(
        "w", newline="", encoding="utf-8"
    ) as ff:
        sum_writer = csv.DictWriter(sf, fieldnames=_SUMMARY_FIELDS, extrasaction="ignore")
        fail_writer = csv.DictWriter(ff, fieldnames=_FAILED_FIELDS, extrasaction="ignore")
        sum_writer.writeheader()
        fail_writer.writeheader()

        for p in xyz_paths:
            name = p.stem
            stem = name.split("_")[-1]
            json_path = json_dir / f"{name}.json"
            if not json_path.is_file():
                fail_writer.writerow({
                    "stem": stem,
                    "error_type": "MissingOutput",
                    "error": "no json/<name>.json (task failed or never ran — check HQ logs)",
                    "elapsed_s": "",
                })
                n_fail += 1
                continue
            try:
                with json_path.open(encoding="utf-8") as fh:
                    data = json.load(fh)
                sum_writer.writerow({
                    "stem": stem,
                    "formula": data["formula"],
                    "tae_Ha": f"{float(data['tae']):.3f}",
                    "Hf_298K_kJ": f"{float(data['Hf_298K_kJ']):.3f}",
                    "LFL_percent": f"{float(data['LFL_percent']):.3f}",
                    "UFL_percent": f"{float(data['UFL_percent']):.3f}",
                    "n_freqs": len(data.get("freqs") or []),
                    "elapsed_s": "",
                })
                n_ok += 1
            except (json.JSONDecodeError, KeyError, ValueError) as exc:
                fail_writer.writerow({
                    "stem": stem,
                    "error_type": type(exc).__name__,
                    "error": " ".join(str(exc).split()) or "malformed JSON",
                    "elapsed_s": "",
                })
                n_fail += 1

    print(
        f"Collected {len(xyz_paths)} molecules from {input_list_path.name}: "
        f"{n_ok} ok, {n_fail} failed."
    )
    print(f"Summary: {summary_csv}")
    if n_fail:
        print(f"Failures: {failed_csv}")
