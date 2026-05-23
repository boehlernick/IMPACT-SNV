#!/usr/bin/env python3
"""Command-line interface for IMPACT-SNV."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

DEFAULT_VERSION = "1.0.0a0"
QC_MODES = ("warn", "strict", "off")


@dataclass(frozen=True)
class FinalizeGdsArgs:
    input_dir: Path
    gene_list: Path
    out_dir: Path
    qc_mode: str
    samples: Optional[list[str]]
    rscript: str
    force: bool
    no_optimize: bool
    keep_intermediate: bool
    manifest_json: Optional[Path]


def _existing_file(path_text: str) -> Path:
    path = Path(path_text)
    if not path.exists() or not path.is_file():
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return path


def _existing_dir(path_text: str) -> Path:
    path = Path(path_text)
    if not path.exists() or not path.is_dir():
        raise argparse.ArgumentTypeError(f"Directory not found: {path}")
    return path


def _path(path_text: str) -> Path:
    return Path(path_text)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="impact-snv",
        description="IMPACT-SNV production CLI for SNV/indel processing through IMPACT-VIS-ready GDS output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"impact-snv {DEFAULT_VERSION}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)
    add_build_gds_parser(subparsers)
    add_finalize_gds_parser(subparsers)
    add_validate_gds_parser(subparsers)
    add_qc_build_parser(subparsers)
    add_placeholder_parsers(subparsers)
    return parser


def add_build_gds_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "build-gds",
        help="Build pre-prioritization per-sample GDS files from FAVOR annotation/genotype parquet outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--annotated-dir", required=True, type=_existing_dir, help="FAVOR annotated output directory")
    p.add_argument("--genotypes-dir", required=True, type=_existing_dir, help="FAVOR genotype output directory containing samples.txt")
    p.add_argument("--gene-list", required=True, type=_existing_file, help="GeneList.txt with symbol and globalScore columns")
    p.add_argument("--out-dir", required=True, type=_path, help="Output build directory")
    p.add_argument("--samples", nargs="+", help="Sample IDs to process; defaults to all samples if omitted")
    p.add_argument("--all-samples", action="store_true", help="Process all samples from genotypes_dir/samples.txt")
    p.add_argument("--chromosomes", nargs="+", default=["1-22", "X", "Y"], help="Chromosomes to process; supports ranges like 1-22")
    p.add_argument("--dosage-threshold", type=float, default=0.0, help="Keep variants with dosage > threshold")
    p.add_argument("--qc-mode", choices=QC_MODES, default="warn", help="QC behavior: warn, strict, or off")
    p.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    p.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    p.add_argument("--python", default=sys.executable, help="Python executable used for flattening subprocesses")
    p.add_argument("--rscript", default="Rscript", help="Rscript executable")
    p.add_argument("--no-optimize", action="store_true", help="Pass --no-optimize to GDS writer")
    p.add_argument("--keep-temp-vcf", action="store_true", help="Keep temporary VCF files from GDS writing")
    p.add_argument("--write-preview-csv", action="store_true", help="Write preview CSVs for flat parquet outputs")
    p.add_argument("--skip-chrom-gds", action="store_true", help="Skip per-chromosome GDS files and only write per-sample merged GDS")
    p.add_argument("--no-merge-per-sample", dest="merge_per_sample", action="store_false", help="Do not create merged per-sample pre-prioritization GDS files")
    p.set_defaults(merge_per_sample=True)
    p.add_argument("--manifest-json", type=_path, help="Optional build manifest path")
    p.set_defaults(func=cmd_build_gds)


def add_finalize_gds_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "finalize-gds",
        help="Score pre-prioritization GDS files and make final IMPACT-VIS-ready *_SNV_IMPACT.gds outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input-dir", required=True, type=_existing_dir, help="Directory containing *_SNV_IMPACT.preprioritization.gds files")
    p.add_argument("--gene-list", required=True, type=_existing_file, help="GeneList.txt with symbol and globalScore columns")
    p.add_argument("--out-dir", required=True, type=_path, help="Directory for final <sample_id>_SNV_IMPACT.gds files")
    p.add_argument("--qc-mode", choices=QC_MODES, default="warn")
    p.add_argument("--samples", nargs="+", default=None, help="Optional sample IDs to finalize")
    p.add_argument("--rscript", default="Rscript")
    p.add_argument("--manifest-json", type=_path, default=None)
    p.add_argument("--force", action="store_true")
    p.add_argument("--no-optimize", action="store_true")
    p.add_argument("--keep-intermediate", action="store_true")
    p.set_defaults(func=cmd_finalize_gds)


def add_validate_gds_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser("validate-gds", help="Validate final IMPACT-VIS-ready *_SNV_IMPACT.gds files", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--gds", type=_existing_file, help="Single final GDS file to validate")
    group.add_argument("--input-dir", type=_existing_dir, help="Directory of final GDS files to validate")
    p.add_argument("--rscript", default="Rscript")
    p.add_argument("--qc-mode", choices=QC_MODES, default="warn")
    p.set_defaults(func=cmd_validate_gds)


def add_qc_build_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser("qc-build", help="Create QC reports from build_manifest.json", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--manifest", required=True, type=_existing_file)
    p.add_argument("--out-dir", required=True, type=_path)
    p.add_argument("--qc-mode", choices=QC_MODES, default="warn")
    p.set_defaults(func=cmd_qc_build)


def add_placeholder_parsers(subparsers: argparse._SubParsersAction) -> None:
    planned = {
        "run": "Run the full workflow from VCF manifest to IMPACT-VIS-ready GDS outputs.",
        "merge": "Merge individual per-sample VCFs into a normalized multi-sample VCF.",
        "favor-ingest": "Run FAVOR CLI ingest.",
        "favor-annotate": "Run FAVOR CLI annotation.",
    }
    for name, help_text in planned.items():
        p = subparsers.add_parser(name, help=f"[planned] {help_text}")
        p.set_defaults(func=cmd_not_implemented, planned_command=name)


def cmd_build_gds(args: argparse.Namespace) -> int:
    try:
        from impact_snv.gds.build import run_build_gds
        from impact_snv.gds.build import package_path, resource_script
    except ModuleNotFoundError as exc:
        print("ERROR: impact_snv.gds.build is not importable", file=sys.stderr)
        raise SystemExit(2) from exc
    # Resolve internal scripts/resources here, after package import succeeds.
    args.flatten_script = package_path("impact_snv.gds.flatten")
    args.gds_writer = resource_script("favor_flat_to_seqarray_gds.R")
    if not args.all_samples and not args.samples:
        args.all_samples = True
    return int(run_build_gds(args) or 0)


def cmd_finalize_gds(args: argparse.Namespace) -> int:
    normalized = FinalizeGdsArgs(
        input_dir=args.input_dir,
        gene_list=args.gene_list,
        out_dir=args.out_dir,
        qc_mode=args.qc_mode,
        samples=args.samples,
        rscript=args.rscript,
        force=args.force,
        no_optimize=args.no_optimize,
        keep_intermediate=args.keep_intermediate,
        manifest_json=args.manifest_json,
    )
    from impact_snv.gds.finalize import run_finalize_gds
    return int(run_finalize_gds(normalized) or 0)


def cmd_validate_gds(args: argparse.Namespace) -> int:
    from impact_snv.gds.validate import run_validate_gds
    return int(run_validate_gds(args) or 0)


def cmd_qc_build(args: argparse.Namespace) -> int:
    from impact_snv.qc.build_qc import run_qc_build
    return int(run_qc_build(args) or 0)


def cmd_not_implemented(args: argparse.Namespace) -> int:
    command = getattr(args, "planned_command", "unknown")
    print(f"impact-snv {command} is planned for IMPACT-SNV v1.0.0 but is not implemented yet.", file=sys.stderr)
    return 2


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
