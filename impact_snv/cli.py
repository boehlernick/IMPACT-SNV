#!/usr/bin/env python3
"""Command-line interface for IMPACT-SNV.

This module defines the production CLI entry point exposed by pyproject.toml:

    impact-snv = "impact_snv.cli:main"

The CLI is intentionally thin. Heavy workflow logic should live in package
modules such as impact_snv.gds.finalize and impact_snv.qc.build_qc so that the
same code can be reused by local execution, tests, and DNAnexus wrappers.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence


DEFAULT_VERSION = "1.0.0a0"
QC_MODES = ("warn", "strict", "off")
DEFAULT_CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]


@dataclass(frozen=True)
class FinalizeGdsArgs:
    """Normalized arguments for the finalize-gds subcommand."""

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


def add_common_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--version",
        action="version",
        version=f"impact-snv {DEFAULT_VERSION}",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="impact-snv",
        description=(
            "IMPACT-SNV production CLI for SNV/indel processing, FAVOR-derived "
            "annotation handling, GDS creation, IMPACT prioritization, and "
            "IMPACT-VIS-ready output generation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_common_options(parser)

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="COMMAND",
        required=True,
    )

    add_finalize_gds_parser(subparsers)
    add_validate_gds_parser(subparsers)
    add_placeholder_parsers(subparsers)

    return parser


def add_finalize_gds_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "finalize-gds",
        help="Score pre-prioritization GDS files and make final IMPACT-VIS-ready *_SNV_IMPACT.gds outputs.",
        description=(
            "Finalize per-sample pre-prioritization GDS files by running IMPACT "
            "prioritization, adding IMPACT-VIS compatibility nodes such as "
            "FunctionalAnnotation/VarInfo, validating final nodes, and writing "
            "a finalize_manifest.json."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--input-dir",
        required=True,
        type=_existing_dir,
        help="Directory containing *_SNV_IMPACT.preprioritization.gds files.",
    )
    p.add_argument(
        "--gene-list",
        required=True,
        type=_existing_file,
        help="GeneList.txt with at least symbol and globalScore columns.",
    )
    p.add_argument(
        "--out-dir",
        required=True,
        type=_path,
        help="Directory where final <sample_id>_SNV_IMPACT.gds files will be written.",
    )
    p.add_argument(
        "--qc-mode",
        choices=QC_MODES,
        default="warn",
        help=(
            "QC behavior. 'warn' records warnings and continues; 'strict' may "
            "raise failures in downstream validation; 'off' performs only minimal checks."
        ),
    )
    p.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help=(
            "Optional sample IDs to finalize. If omitted, all "
            "*_SNV_IMPACT.preprioritization.gds files in --input-dir are processed."
        ),
    )
    p.add_argument(
        "--rscript",
        default="Rscript",
        help="Rscript executable to use for R resource scripts.",
    )
    p.add_argument(
        "--manifest-json",
        type=_path,
        default=None,
        help="Optional manifest path. Defaults to <out-dir>/finalize_manifest.json.",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing final output GDS files.",
    )
    p.add_argument(
        "--no-optimize",
        action="store_true",
        help="Pass --no-optimize to R finalization resources where supported.",
    )
    p.add_argument(
        "--keep-intermediate",
        action="store_true",
        help=(
            "Keep intermediate scored GDS files produced between prioritization "
            "and IMPACT-VIS compatibility injection."
        ),
    )
    p.set_defaults(func=cmd_finalize_gds)


def add_validate_gds_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "validate-gds",
        help="Validate final IMPACT-VIS-ready *_SNV_IMPACT.gds files.",
        description=(
            "Validate final GDS files for required IMPACT-SNV/IMPACT-VIS nodes. "
            "This command will be backed by impact_snv.gds.validate."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--gds",
        type=_existing_file,
        help="Single final *_SNV_IMPACT.gds file to validate.",
    )
    group.add_argument(
        "--input-dir",
        type=_existing_dir,
        help="Directory of final *_SNV_IMPACT.gds files to validate.",
    )
    p.add_argument(
        "--rscript",
        default="Rscript",
        help="Rscript executable, if validation requires an R backend.",
    )
    p.add_argument(
        "--qc-mode",
        choices=QC_MODES,
        default="warn",
        help="QC behavior for validation warnings.",
    )
    p.set_defaults(func=cmd_validate_gds)


def add_placeholder_parsers(subparsers: argparse._SubParsersAction) -> None:
    """Register planned v1.0.0 commands with helpful not-yet-implemented messages.

    These placeholders make the intended production CLI visible early while we
    implement one stable command at a time.
    """
    planned = {
        "run": "Run the full local IMPACT-SNV workflow from VCF manifest to IMPACT-VIS-ready GDS outputs.",
        "merge": "Merge individual per-sample VCFs into a normalized multi-sample VCF.",
        "favor-ingest": "Run FAVOR CLI ingest for merged/normalized SNV/indel VCF inputs.",
        "favor-annotate": "Run FAVOR CLI annotation for ingested variant data.",
        "build-gds": "Build pre-prioritization per-sample GDS files from FAVOR annotation/genotype parquet outputs.",
        "qc-build": "Summarize build manifests and sample/chromosome row-count warnings.",
    }
    for name, help_text in planned.items():
        p = subparsers.add_parser(
            name,
            help=f"[planned] {help_text}",
            description=f"Planned v1.0.0 command: {help_text}",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        p.set_defaults(func=cmd_not_implemented, planned_command=name)


def cmd_finalize_gds(args: argparse.Namespace) -> int:
    """Dispatch finalize-gds to impact_snv.gds.finalize.

    The import is intentionally lazy so that `impact-snv --help` works even if
    optional runtime dependencies for finalization are not yet installed.
    """
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

    try:
        from impact_snv.gds.finalize import run_finalize_gds
    except ModuleNotFoundError as exc:
        print(
            "ERROR: impact_snv.gds.finalize is not implemented or not importable yet.\n"
            "Next step: create impact_snv/gds/finalize.py with run_finalize_gds(args).",
            file=sys.stderr,
        )
        raise SystemExit(2) from exc

    return int(run_finalize_gds(normalized) or 0)


def cmd_validate_gds(args: argparse.Namespace) -> int:
    try:
        from impact_snv.gds.validate import run_validate_gds
    except ModuleNotFoundError as exc:
        print(
            "ERROR: impact_snv.gds.validate is not implemented or not importable yet.\n"
            "Next step: create impact_snv/gds/validate.py with run_validate_gds(args).",
            file=sys.stderr,
        )
        raise SystemExit(2) from exc

    return int(run_validate_gds(args) or 0)


def cmd_not_implemented(args: argparse.Namespace) -> int:
    command = getattr(args, "planned_command", "unknown")
    print(
        f"impact-snv {command!s} is planned for IMPACT-SNV v1.0.0 but is not implemented yet.\n"
        "Current implementation focus: finalize-gds, then build-gds, then full run orchestration.",
        file=sys.stderr,
    )
    return 2


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
