from __future__ import annotations

import argparse

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
        help="In batch mode, skip inputs whose output JSON already exists.",
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

    if args.xyz is None and args.input_list is None:
        parser.error("either an xyz positional argument or --input-list is required.")
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

    run_pipeline(cfg)


if __name__ == "__main__":
    main()
