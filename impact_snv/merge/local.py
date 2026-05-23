#!/usr/bin/env python3
"""Local VCF merge CLI backend for IMPACT-SNV.

This module is used both as a standalone module and by impact_snv.cli for the
production command:

    impact-snv merge
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Sequence

from impact_snv.merge.core import MergeConfig, MergeResult, VCFMergeError, merge_result_manifest, merge_vcf_samples

QC_MODES = {"warn", "strict", "off"}
NORMALIZATION_MODES = {"auto", "always", "never"}
VERSION = "1.0.0a0"


@dataclass(frozen=True)
class MergeCliArgs:
    vcfs: list[Path]
    out_vcf: Path
    reference_fasta: Optional[Path]
    reference_build: str
    normalization_mode: str
    preflight_records: int
    threads: int
    work_dir: Optional[Path]
    keep_intermediates: bool
    force: bool
    force_samples: bool
    bcftools: Optional[Path]
    bgzip: Optional[Path]
    extra_merge_arg: list[str]
    command_log: Optional[Path]
    manifest_json: Optional[Path]
    qc_mode: str
    quiet: bool


def positive_int(v: str) -> int:
    try:
        x = int(v)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Expected integer, got {v!r}") from e
    if x < 1:
        raise argparse.ArgumentTypeError("Value must be >= 1")
    return x


def load_vcfs_from_manifest(path: Path) -> list[Path]:
    """Load VCF paths from a TSV manifest.

    Accepted columns:
      - vcf_path
      - vcf
      - path

    A sample_id column is allowed but not currently used for reheadering.
    """
    if not path.exists():
        raise FileNotFoundError(f"VCF manifest not found: {path}")
    text = path.read_text(encoding="utf-8").splitlines()
    rows = [line for line in text if line.strip() and not line.lstrip().startswith("#")]
    if not rows:
        raise ValueError(f"VCF manifest is empty: {path}")
    header = rows[0].rstrip("\n").split("\t")
    candidates = ["vcf_path", "vcf", "path"]
    idx = next((header.index(c) for c in candidates if c in header), None)
    if idx is None:
        raise ValueError(f"VCF manifest must contain one of columns {candidates}; found {header}")
    vcfs = [Path(row.split("\t")[idx]) for row in rows[1:] if len(row.split("\t")) > idx]
    if len(vcfs) < 2:
        raise ValueError("VCF manifest must contain at least two VCF paths")
    return vcfs


def make_config(args: MergeCliArgs) -> MergeConfig:
    return MergeConfig(
        reference_fasta=args.reference_fasta,
        normalization_mode=args.normalization_mode,
        preflight_records=args.preflight_records,
        threads=args.threads,
        keep_intermediates=args.keep_intermediates,
        force_samples=args.force_samples,
        overwrite=args.force,
        work_dir=args.work_dir,
        bcftools_path=args.bcftools,
        bgzip_path=args.bgzip,
        extra_merge_args=tuple(args.extra_merge_arg or []),
    )


def write_text(path: Path, text: str) -> None:
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_merge_manifest(path: Path, result: MergeResult, args: MergeCliArgs, status: str = "ok", message: str = "") -> None:
    payload = {
        "command": "impact-snv merge",
        "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "message": message,
        "reference_build": args.reference_build,
        "threads": args.threads,
        "force_samples": args.force_samples,
        "extra_merge_args": list(args.extra_merge_arg or []),
        **merge_result_manifest(result, qc_mode=args.qc_mode, reference_build=args.reference_build),
    }
    write_text(path, json.dumps(payload, indent=2, sort_keys=True))


def run_merge(args: MergeCliArgs) -> int:
    if args.reference_build != "GRCh38":
        print("ERROR: IMPACT-SNV v1.0.0 supports only --reference-build GRCh38", file=sys.stderr)
        return 2
    if len(args.vcfs) < 2:
        print("ERROR: merge requires at least two input VCF files", file=sys.stderr)
        return 2
    if args.normalization_mode == "always" and args.reference_fasta is None:
        print("ERROR: --normalization-mode always requires --reference-fasta", file=sys.stderr)
        return 2

    try:
        result = merge_vcf_samples(args.vcfs, args.out_vcf, make_config(args))
    except VCFMergeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("ERROR: interrupted by user", file=sys.stderr)
        return 130

    manifest_path = args.manifest_json or Path(str(result.merged_vcf) + ".merge_manifest.json")
    write_merge_manifest(manifest_path, result, args)
    if args.command_log:
        write_text(args.command_log, "\n".join(result.command_log) + "\n")

    if not args.quiet:
        pf = result.normalization_preflight
        print("VCF merge completed successfully")
        print(f"  merged_vcf: {result.merged_vcf}")
        print(f"  merged_index: {result.merged_index or 'not found'}")
        print(f"  samples: {len(result.sample_ids)}")
        print(f"  normalization_mode: {pf.mode}")
        print(f"  normalization_performed: {pf.should_normalize}")
        print(f"  manifest: {manifest_path}")
        if pf.reasons:
            print("  normalization_reasons:")
            for reason in pf.reasons:
                print(f"   - {reason}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Merge VCF/VCF.GZ/BCF files into an indexed multi-sample VCF.GZ")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--vcfs", nargs="+", type=Path, help="Input VCF/VCF.GZ/BCF files; at least two")
    src.add_argument("--vcf-manifest", type=Path, help="TSV manifest with a vcf_path/vcf/path column")
    p.add_argument("--out-vcf", required=True, type=Path, help="Output merged VCF path")
    p.add_argument("--reference-fasta", type=Path, help="Reference FASTA for bcftools norm; requires .fai")
    p.add_argument("--reference-build", default="GRCh38", choices=["GRCh38"])
    p.add_argument("--normalization-mode", choices=sorted(NORMALIZATION_MODES), default="auto")
    p.add_argument("--preflight-records", type=positive_int, default=10000)
    p.add_argument("--threads", type=positive_int, default=1)
    p.add_argument("--work-dir", type=Path)
    p.add_argument("--keep-intermediates", action="store_true")
    p.add_argument("--force", action="store_true", help="Overwrite existing output/intermediates")
    p.add_argument("--force-samples", action="store_true", help="Allow duplicate sample IDs via bcftools --force-samples")
    p.add_argument("--bcftools", type=Path)
    p.add_argument("--bgzip", type=Path)
    p.add_argument("--extra-merge-arg", action="append", default=[])
    p.add_argument("--command-log", type=Path)
    p.add_argument("--manifest-json", type=Path)
    p.add_argument("--qc-mode", choices=sorted(QC_MODES), default="warn")
    p.add_argument("--quiet", action="store_true")
    return p


def args_from_namespace(ns: argparse.Namespace) -> MergeCliArgs:
    vcfs = [Path(x) for x in ns.vcfs] if ns.vcfs else load_vcfs_from_manifest(Path(ns.vcf_manifest))
    return MergeCliArgs(
        vcfs=vcfs,
        out_vcf=Path(ns.out_vcf),
        reference_fasta=Path(ns.reference_fasta) if ns.reference_fasta else None,
        reference_build=ns.reference_build,
        normalization_mode=ns.normalization_mode,
        preflight_records=ns.preflight_records,
        threads=ns.threads,
        work_dir=Path(ns.work_dir) if ns.work_dir else None,
        keep_intermediates=ns.keep_intermediates,
        force=ns.force,
        force_samples=ns.force_samples,
        bcftools=Path(ns.bcftools) if ns.bcftools else None,
        bgzip=Path(ns.bgzip) if ns.bgzip else None,
        extra_merge_arg=list(ns.extra_merge_arg or []),
        command_log=Path(ns.command_log) if ns.command_log else None,
        manifest_json=Path(ns.manifest_json) if ns.manifest_json else None,
        qc_mode=ns.qc_mode,
        quiet=ns.quiet,
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    ns = build_parser().parse_args(argv)
    return run_merge(args_from_namespace(ns))


if __name__ == "__main__":
    raise SystemExit(main())
