from __future__ import annotations

import argparse
import json
import time
import traceback
from pathlib import Path

from .config import build_config, validate_config
from .stages import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the flammability pipeline. Default mode 'ml' computes TAE and "
        "frequencies from an XYZ file using skala+hip. Mode 'ab' reads them from a small YAML."
    )
    parser.add_argument(
        "xyz",
        nargs="?",
        default=None,
        help="Path to a single XYZ geometry file. Required unless --input-list is set.",
    )
    parser.add_argument(
        "--input-list",
        default=None,
        help="Text file with one XYZ path per line; runs the pipeline once per molecule "
        "in batch mode. Mutually exclusive with the positional xyz argument.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel workers in batch mode (default: 1). Only used with --input-list.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip molecules whose output JSON already exists. Works in both single-molecule "
        "and batch mode; useful as a HyperQueue retry guard.",
    )
    parser.add_argument(
        "--collect-summary",
        metavar="XYZ_LIST",
        default=None,
        help="Post-processing: scan <output_dir>/json/ and rebuild SUMMARY.csv from the "
        "per-molecule JSONs. Molecules listed in XYZ_LIST but missing a JSON are written "
        "to FAILED.csv. Use this after a HyperQueue per-molecule batch finishes.",
    )
    parser.add_argument(
        "--mode",
        choices=["ml", "ab"],
        default="ml",
        help="ml: compute TAE/freqs from XYZ via skala+hip (default). "
        "ab: read TAE/freqs from --yaml.",
    )
    parser.add_argument(
        "--yaml",
        dest="case_yaml",
        default=None,
        help="Case-data YAML with 'tae' (Hartree) and 'freqs' (cm^-1); required for --mode ab.",
    )
    parser.add_argument("--output-dir", default=None, help="Override output directory.")
    parser.add_argument("--no-plot", action="store_true", help="Skip writing the phase-diagram PDF.")
    parser.add_argument(
        "--hip-checkpoint",
        default=None,
        help="Path to hip checkpoint (or set HIP_CKPT). "
        "If unset, hip_v3_cf.ckpt is downloaded from huggingface.co/andreasburger/hip on first use.",
    )
    parser.add_argument(
        "--hip-device",
        default="auto",
        choices=["auto", "cpu", "cuda"],
        help="hip torch device. 'auto' picks cuda for molecules >10 heavy atoms, "
        "cpu otherwise (default: auto).",
    )
    parser.add_argument(
        "--skala-device",
        default="auto",
        choices=["auto", "cpu", "cuda"],
        help="skala torch device. 'auto' picks cuda for molecules >10 heavy atoms, "
        "cpu otherwise (default: auto).",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate inputs and exit without running the pipeline.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.collect_summary is not None:
        from .batch import collect_summary

        collect_summary(
            input_list=args.collect_summary,
            output_dir=args.output_dir,
        )
        return

    if args.xyz is None and args.input_list is None:
        parser.error(
            "either an xyz positional argument, --input-list, or --collect-summary is required."
        )
    if args.xyz is not None and args.input_list is not None:
        parser.error("--input-list cannot be combined with a positional xyz argument.")

    if args.input_list is not None:
        from .batch import run_batch

        run_batch(
            input_list=args.input_list,
            mode=args.mode,
            case_yaml=args.case_yaml,
            output_dir=args.output_dir,
            plot_map=not args.no_plot,
            hip_checkpoint=args.hip_checkpoint,
            hip_device=args.hip_device,
            skala_device=args.skala_device,
            jobs=args.jobs,
            skip_existing=args.skip_existing,
            validate_only=args.validate_only,
        )
        return

    cfg = build_config(
        xyz_geom=args.xyz,
        mode=args.mode,
        case_yaml=args.case_yaml,
        output_dir=args.output_dir,
        plot_map=not args.no_plot,
        hip_checkpoint=args.hip_checkpoint,
        hip_device=args.hip_device,
        skala_device=args.skala_device,
    )

    errors = validate_config(cfg)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        raise SystemExit(1)

    if args.validate_only:
        print("Inputs are valid.")
        return

    stem = Path(args.xyz).stem
    json_dir = cfg.output_dir / "json"
    ok_json = json_dir / f"{stem}.json"
    failed_json = json_dir / f"{stem}.failed.json"

    # --skip-existing in single mode is the HyperQueue retry guard: any prior
    # verdict (success OR recorded failure) counts as "already attempted" and
    # exits 0 so HQ marks the task done. To redo failures, `rm json/*.failed.json`
    # before resubmitting.
    if args.skip_existing:
        for p in (ok_json, failed_json):
            if p.is_file():
                print(f"--skip-existing: {p} already exists, skipping.")
                return

    # Wrap the pipeline so HQ-mode failures leave a structured record next to
    # the success JSONs. Without this, an exception in run_pipeline only ends
    # up in HQ's per-task stderr — and `hq_logs/job-NNN.stderr` has no naming
    # link back to the xyz stem, making collect_summary's FAILED.csv useless.
    json_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    try:
        run_pipeline(cfg, case_name_override=stem)
    except Exception as exc:
        failed_json.write_text(
            json.dumps(
                {
                    "status": "fail",
                    "error_type": type(exc).__name__,
                    "error": " ".join(str(exc).split()),
                    "traceback": traceback.format_exc(),
                    "xyz_geom": str(cfg.xyz_geom),
                    "elapsed_s": round(time.perf_counter() - t0, 3),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        # Re-raise so the traceback still hits stderr and HQ records the task
        # as failed (a silent exit 0 would make HQ think it succeeded).
        raise


if __name__ == "__main__":
    main()
