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


def positive_int(v: str) -> int:
    try:
        x = int(v)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Expected integer, got {v!r}") from e
    if x < 1:
        raise argparse.ArgumentTypeError("Value must be >= 1")
    return x


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="impact-snv",
        description="IMPACT-SNV production CLI for SNV/indel processing through IMPACT-VIS-ready GDS output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"impact-snv {DEFAULT_VERSION}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)
    add_merge_parser(subparsers)
    add_favor_ingest_parser(subparsers)
    add_favor_annotate_parser(subparsers)
    add_build_gds_parser(subparsers)
    add_finalize_gds_parser(subparsers)
    add_validate_gds_parser(subparsers)
    add_qc_build_parser(subparsers)
    add_placeholder_parsers(subparsers)
    return parser


def add_merge_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser("merge", help="Merge individual per-sample VCFs into an indexed multi-sample VCF.GZ.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--vcfs", nargs="+", type=Path)
    src.add_argument("--vcf-manifest", type=_existing_file)
    p.add_argument("--out-vcf", required=True, type=_path)
    p.add_argument("--reference-fasta", type=_existing_file)
    p.add_argument("--reference-build", default="GRCh38", choices=["GRCh38"])
    p.add_argument("--normalization-mode", choices=["auto", "always", "never"], default="auto")
    p.add_argument("--preflight-records", type=positive_int, default=10000)
    p.add_argument("--threads", type=positive_int, default=1)
    p.add_argument("--work-dir", type=_path)
    p.add_argument("--keep-intermediates", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--force-samples", action="store_true")
    p.add_argument("--bcftools", type=_path)
    p.add_argument("--bgzip", type=_path)
    p.add_argument("--extra-merge-arg", action="append", default=[])
    p.add_argument("--command-log", type=_path)
    p.add_argument("--manifest-json", type=_path)
    p.add_argument("--qc-mode", choices=QC_MODES, default="warn")
    p.add_argument("--quiet", action="store_true")
    p.set_defaults(func=cmd_merge)


def add_favor_ingest_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser("favor-ingest", help="Run FAVOR CLI ingest using the tested IMPACT-SNV contract.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--input-vcf", required=True, type=_existing_file)
    p.add_argument("--out-dir", required=True, type=_path)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--reference-build", default="GRCh38", choices=["GRCh38"])
    p.add_argument("--favor-bin", default="favor")
    p.add_argument("--force", action="store_true")
    p.add_argument("--manifest-json", type=_path)
    p.add_argument("--progress-interval-seconds", type=positive_int, default=60)
    p.add_argument("--quiet", action="store_true")
    p.add_argument("--tail-log-lines", type=int, default=0)
    p.add_argument("--max-progress-log-line-chars", type=positive_int, default=300)
    p.add_argument("--progress-mode", choices=["compact", "normal", "verbose"], default="normal")
    p.set_defaults(func=cmd_favor_ingest)


def add_favor_annotate_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser("favor-annotate", help="Run FAVOR CLI annotate using the tested IMPACT-SNV contract.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--ingested-dir", required=True, type=_existing_dir)
    p.add_argument("--out-dir", required=True, type=_path)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--reference-build", default="GRCh38", choices=["GRCh38"])
    p.add_argument("--favor-bin", default="favor")
    p.add_argument("--force", action="store_true")
    p.add_argument("--manifest-json", type=_path)
    p.add_argument("--progress-interval-seconds", type=positive_int, default=60)
    p.add_argument("--quiet", action="store_true")
    p.add_argument("--tail-log-lines", type=int, default=0)
    p.add_argument("--max-progress-log-line-chars", type=positive_int, default=300)
    p.add_argument("--progress-mode", choices=["compact", "normal", "verbose"], default="normal")
    p.set_defaults(func=cmd_favor_annotate)


def add_build_gds_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser("build-gds", help="Build pre-prioritization per-sample GDS files from FAVOR parquet outputs.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--annotated-dir", required=True, type=_existing_dir)
    p.add_argument("--genotypes-dir", required=True, type=_existing_dir)
    p.add_argument("--gene-list", required=True, type=_existing_file)
    p.add_argument("--out-dir", required=True, type=_path)
    p.add_argument("--samples", nargs="+")
    p.add_argument("--all-samples", action="store_true")
    p.add_argument("--chromosomes", nargs="+", default=["1-22", "X", "Y"])
    p.add_argument("--dosage-threshold", type=float, default=0.0)
    p.add_argument("--qc-mode", choices=QC_MODES, default="warn")
    p.add_argument("--force", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--python", default=sys.executable)
    p.add_argument("--rscript", default="Rscript")
    p.add_argument("--no-optimize", action="store_true")
    p.add_argument("--keep-temp-vcf", action="store_true")
    p.add_argument("--write-preview-csv", action="store_true")
    p.add_argument("--skip-chrom-gds", action="store_true")
    p.add_argument("--no-merge-per-sample", dest="merge_per_sample", action="store_false")
    p.set_defaults(merge_per_sample=True)
    p.add_argument("--manifest-json", type=_path)
    p.set_defaults(func=cmd_build_gds)


def add_finalize_gds_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser("finalize-gds", help="Score pre-prioritization GDS files and make final IMPACT-VIS-ready outputs.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--input-dir", required=True, type=_existing_dir)
    p.add_argument("--gene-list", required=True, type=_existing_file)
    p.add_argument("--out-dir", required=True, type=_path)
    p.add_argument("--qc-mode", choices=QC_MODES, default="warn")
    p.add_argument("--samples", nargs="+", default=None)
    p.add_argument("--rscript", default="Rscript")
    p.add_argument("--manifest-json", type=_path, default=None)
    p.add_argument("--force", action="store_true")
    p.add_argument("--no-optimize", action="store_true")
    p.add_argument("--keep-intermediate", action="store_true")
    p.set_defaults(func=cmd_finalize_gds)


def add_validate_gds_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser("validate-gds", help="Validate final IMPACT-VIS-ready GDS files", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--gds", type=_existing_file)
    group.add_argument("--input-dir", type=_existing_dir)
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
    p = subparsers.add_parser("run", help="[planned] Run the full workflow from VCF manifest to IMPACT-VIS-ready GDS outputs.")
    p.set_defaults(func=cmd_not_implemented, planned_command="run")


def cmd_merge(args: argparse.Namespace) -> int:
    from impact_snv.merge.local import MergeCliArgs, load_vcfs_from_manifest, run_merge
    vcfs = [Path(x) for x in args.vcfs] if args.vcfs else load_vcfs_from_manifest(Path(args.vcf_manifest))
    return run_merge(MergeCliArgs(vcfs, args.out_vcf, args.reference_fasta, args.reference_build, args.normalization_mode, args.preflight_records, args.threads, args.work_dir, args.keep_intermediates, args.force, args.force_samples, args.bcftools, args.bgzip, list(args.extra_merge_arg or []), args.command_log, args.manifest_json, args.qc_mode, args.quiet))


def cmd_favor_ingest(args: argparse.Namespace) -> int:
    from impact_snv.favor.ingest import run_favor_ingest
    return int(run_favor_ingest(args) or 0)


def cmd_favor_annotate(args: argparse.Namespace) -> int:
    from impact_snv.favor.annotate import run_favor_annotate
    return int(run_favor_annotate(args) or 0)


def cmd_build_gds(args: argparse.Namespace) -> int:
    from impact_snv.gds.build import package_path, resource_script, run_build_gds
    args.flatten_script = package_path("impact_snv.gds.flatten")
    args.gds_writer = resource_script("favor_flat_to_seqarray_gds.R")
    if not args.all_samples and not args.samples:
        args.all_samples = True
    return int(run_build_gds(args) or 0)


def cmd_finalize_gds(args: argparse.Namespace) -> int:
    normalized = FinalizeGdsArgs(args.input_dir, args.gene_list, args.out_dir, args.qc_mode, args.samples, args.rscript, args.force, args.no_optimize, args.keep_intermediate, args.manifest_json)
    from impact_snv.gds.finalize import run_finalize_gds
    return int(run_finalize_gds(normalized) or 0)


def cmd_validate_gds(args: argparse.Namespace) -> int:
    from impact_snv.gds.validate import run_validate_gds
    return int(run_validate_gds(args) or 0)


def cmd_qc_build(args: argparse.Namespace) -> int:
    from impact_snv.qc.build_qc import run_qc_build
    return int(run_qc_build(args) or 0)


def cmd_not_implemented(args: argparse.Namespace) -> int:
    print(f"impact-snv {getattr(args, 'planned_command', 'unknown')} is planned but not implemented yet.", file=sys.stderr)
    return 2


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
