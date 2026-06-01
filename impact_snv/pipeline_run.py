#!/usr/bin/env python3
"""Local end-to-end orchestration for the supported IMPACT-SNV package path."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Sequence

from impact_snv import __version__ as VERSION
from impact_snv.favor.annotate import BACKEND_FAVOR_CLI, run_favor_annotate
from impact_snv.favor.common import validate_reference_build
from impact_snv.favor.ingest import run_favor_ingest
from impact_snv.gds.build import package_path, resource_script, run_build_gds
from impact_snv.gds.finalize import run_finalize_gds
from impact_snv.gds.validate import validate_gds_dir, write_validation_manifest
from impact_snv.gene_lists import read_samples_manifest, run_build_gene_lists
from impact_snv.genotypes.extract import run_extract_genotypes
from impact_snv.merge.core import (
    VCFMergeError,
    assert_unique_sample_ids,
    get_sample_ids,
    resolve_executable,
)
from impact_snv.merge.local import MergeCliArgs, load_vcfs_from_manifest, run_merge
from impact_snv.merge.sanitize import run_sanitize_vcfs
from impact_snv.qc.build_qc import run_qc_build


@dataclass(frozen=True)
class PipelineStepRecord:
    name: str
    status: str
    message: str = ""
    manifest: Optional[str] = None
    outputs: dict[str, str] = field(default_factory=dict)


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


class PipelineRunError(RuntimeError):
    """Raised when a pipeline step fails."""


def issue(message: str) -> str:
    return message.strip()


def write_run_manifest(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def collect_input_vcfs(samples_manifest: Path) -> tuple[list[str], list[Path]]:
    records = read_samples_manifest(samples_manifest)
    vcfs: list[Path] = []
    for record in records:
        vcf_path = Path(record.vcf_path).expanduser().resolve()
        if not record.vcf_path:
            raise ValueError(f"samples manifest row for {record.sample_id!r} is missing vcf_path")
        if not vcf_path.exists() or not vcf_path.is_file():
            raise FileNotFoundError(f"VCF path for sample {record.sample_id!r} does not exist: {vcf_path}")
        vcfs.append(vcf_path)
    if len(vcfs) < 2:
        raise ValueError("impact-snv run requires at least two input VCF files")
    return [record.sample_id for record in records], vcfs


def detect_duplicate_input_sample_ids(vcfs: Sequence[Path], *, bcftools_path: Optional[Path]) -> tuple[bool, dict[str, list[str]]]:
    command_log: list[str] = []
    bcftools = resolve_executable("bcftools", bcftools_path)
    sample_ids_by_file = {vcf: tuple(get_sample_ids(vcf, bcftools=bcftools, command_log=command_log)) for vcf in vcfs}
    duplicates: dict[str, list[str]] = {}
    try:
        assert_unique_sample_ids(sample_ids_by_file)
        return False, duplicates
    except VCFMergeError:
        for vcf, sample_ids in sample_ids_by_file.items():
            for sample_id in sample_ids:
                duplicates.setdefault(sample_id, []).append(str(vcf))
        duplicates = {sample_id: paths for sample_id, paths in duplicates.items() if len(paths) > 1}
        return True, duplicates


def resolve_gene_source(args: Any, gene_lists_dir: Path) -> tuple[Optional[Path], Optional[Path], Optional[Path]]:
    if getattr(args, "phenotypes_manifest", None):
        manifest_tsv = gene_lists_dir / "sample_gene_lists.tsv"
        manifest_json = gene_lists_dir / "sample_gene_lists.json"
        gene_args = argparse.Namespace(
            samples_manifest=args.samples_manifest,
            phenotypes_manifest=args.phenotypes_manifest,
            out_dir=gene_lists_dir,
            api_url=args.api_url,
            max_search_hits=args.max_search_hits,
            page_size=args.page_size,
            manifest_tsv=manifest_tsv,
            manifest_json=manifest_json,
            force=args.force,
            quiet=args.quiet,
        )
        rc = run_build_gene_lists(gene_args)
        if rc != 0:
            raise PipelineRunError("build-gene-lists failed")
        return None, manifest_tsv, manifest_json
    gene_list = Path(args.gene_list).expanduser().resolve()
    return gene_list, None, None


def run_pipeline(args: Any) -> int:
    reference_build = getattr(args, "reference_build", "GRCh38")
    validate_reference_build(reference_build)

    samples_manifest = Path(args.samples_manifest).expanduser().resolve()
    phenotypes_manifest = (
        Path(args.phenotypes_manifest).expanduser().resolve()
        if getattr(args, "phenotypes_manifest", None)
        else None
    )
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(getattr(args, "manifest_json", None) or (out_dir / "run_manifest.json")).expanduser().resolve()

    steps: list[PipelineStepRecord] = []
    status = "ok"
    message = ""

    try:
        case_ids, vcfs = collect_input_vcfs(samples_manifest)
        sanitize_mode = getattr(args, "sanitize_inputs", "auto") or "auto"
        has_duplicates, duplicate_map = detect_duplicate_input_sample_ids(
            vcfs,
            bcftools_path=Path(args.bcftools) if getattr(args, "bcftools", None) else None,
        )

        use_sanitize = sanitize_mode == "always" or (sanitize_mode == "auto" and has_duplicates)
        if sanitize_mode == "never" and has_duplicates:
            raise PipelineRunError(
                issue(
                    "Duplicate sample IDs were detected in the input VCFs while --sanitize-inputs=never. "
                    "Use --sanitize-inputs auto or always, or sanitize the files before running impact-snv run."
                )
            )

        if use_sanitize:
            sanitize_dir = out_dir / "sanitize"
            sanitize_manifest_tsv = sanitize_dir / "sanitized_vcfs.tsv"
            sanitize_manifest_json = sanitize_dir / "sanitize_manifest.json"
            sanitize_args = argparse.Namespace(
                vcfs=vcfs,
                out_dir=sanitize_dir,
                prefix_regex=args.sanitize_prefix_regex,
                threads=args.threads,
                force=args.force,
                bcftools=args.bcftools,
                bgzip=args.bgzip,
                manifest_tsv=sanitize_manifest_tsv,
                manifest_json=sanitize_manifest_json,
                command_log=sanitize_dir / "sanitize_commands.log",
                quiet=args.quiet,
            )
            rc = run_sanitize_vcfs(sanitize_args)
            if rc != 0:
                raise PipelineRunError("sanitize-vcfs failed")
            merge_inputs = load_vcfs_from_manifest(sanitize_manifest_tsv)
            steps.append(
                PipelineStepRecord(
                    name="sanitize-vcfs",
                    status="ok",
                    manifest=str(sanitize_manifest_json),
                    outputs={"manifest_tsv": str(sanitize_manifest_tsv), "out_dir": str(sanitize_dir)},
                )
            )
        else:
            merge_inputs = vcfs
            steps.append(
                PipelineStepRecord(
                    name="sanitize-vcfs",
                    status="skipped",
                    message="Input sample IDs were already unique; sanitization was not required.",
                )
            )

        merge_dir = out_dir / "merge"
        merge_dir.mkdir(parents=True, exist_ok=True)
        merged_vcf = merge_dir / f"{args.run_prefix}.merged.vcf.gz"
        merge_manifest = merge_dir / f"{args.run_prefix}.merge_manifest.json"
        merge_args = MergeCliArgs(
            vcfs=list(merge_inputs),
            out_vcf=merged_vcf,
            reference_fasta=Path(args.reference_fasta).expanduser().resolve() if getattr(args, "reference_fasta", None) else None,
            reference_build=reference_build,
            normalization_mode=args.normalization_mode,
            preflight_records=args.preflight_records,
            threads=args.threads,
            work_dir=None,
            keep_intermediates=False,
            force=args.force,
            force_samples=False,
            bcftools=Path(args.bcftools).expanduser().resolve() if getattr(args, "bcftools", None) else None,
            bgzip=Path(args.bgzip).expanduser().resolve() if getattr(args, "bgzip", None) else None,
            extra_merge_arg=list(args.extra_merge_arg or []),
            command_log=merge_dir / f"{args.run_prefix}.merge_commands.log",
            manifest_json=merge_manifest,
            qc_mode=args.qc_mode,
            quiet=args.quiet,
        )
        rc = run_merge(merge_args)
        if rc != 0:
            raise PipelineRunError("merge failed")
        steps.append(
            PipelineStepRecord(
                name="merge",
                status="ok",
                manifest=str(merge_manifest),
                outputs={"merged_vcf": str(merged_vcf)},
            )
        )

        favor_dir = out_dir / "favor"
        ingest_manifest = favor_dir / f"{args.run_prefix}.ingest_manifest.json"
        ingest_args = argparse.Namespace(
            input_vcf=merged_vcf,
            out_dir=favor_dir,
            out_prefix=args.run_prefix,
            reference_build=reference_build,
            favor_bin=args.favor_bin,
            force=args.force,
            manifest_json=ingest_manifest,
            progress_interval_seconds=args.progress_interval_seconds,
            quiet=args.quiet,
            tail_log_lines=args.tail_log_lines,
            max_progress_log_line_chars=args.max_progress_log_line_chars,
            progress_mode=args.progress_mode,
        )
        rc = run_favor_ingest(ingest_args)
        if rc != 0:
            raise PipelineRunError("favor-ingest failed")
        ingested_dir = favor_dir / f"{args.run_prefix}.ingested"
        steps.append(
            PipelineStepRecord(
                name="favor-ingest",
                status="ok",
                manifest=str(ingest_manifest),
                outputs={"ingested_dir": str(ingested_dir)},
            )
        )

        annotate_manifest = favor_dir / f"{args.run_prefix}.annotate_manifest.json"
        annotate_args = argparse.Namespace(
            backend=BACKEND_FAVOR_CLI,
            ingested_dir=ingested_dir,
            out_dir=favor_dir,
            out_prefix=args.run_prefix,
            reference_build=reference_build,
            favor_bin=args.favor_bin,
            force=args.force,
            manifest_json=annotate_manifest,
            progress_interval_seconds=args.progress_interval_seconds,
            quiet=args.quiet,
            tail_log_lines=args.tail_log_lines,
            max_progress_log_line_chars=args.max_progress_log_line_chars,
            progress_mode=args.progress_mode,
            legacy_annotated_dir=None,
            legacy_genotypes_dir=None,
            legacy_stage_mode="symlink",
            skeleton_input_file=None,
            skeleton_input_type="gds",
            skeleton_dry_run=True,
            favor_database_file=None,
            favor_database_version=None,
            skeleton_threads=4,
            skeleton_memory_budget_gb=8,
        )
        rc = run_favor_annotate(annotate_args)
        if rc != 0:
            raise PipelineRunError("favor-annotate failed")
        annotated_dir = favor_dir / f"{args.run_prefix}.annotated"
        steps.append(
            PipelineStepRecord(
                name="favor-annotate",
                status="ok",
                manifest=str(annotate_manifest),
                outputs={"annotated_dir": str(annotated_dir)},
            )
        )

        genotypes_dir = out_dir / "genotypes"
        genotype_manifest = genotypes_dir / "genotype_manifest.json"
        extract_args = argparse.Namespace(
            input_vcf=merged_vcf,
            out_dir=genotypes_dir,
            reference_build=reference_build,
            chromosomes=args.chromosomes,
            bcftools=args.bcftools,
            force=args.force,
            qc_mode=args.qc_mode,
            progress_interval_records=args.progress_interval_records,
            manifest_json=genotype_manifest,
        )
        rc = run_extract_genotypes(extract_args)
        if rc != 0:
            raise PipelineRunError("extract-genotypes failed")
        steps.append(
            PipelineStepRecord(
                name="extract-genotypes",
                status="ok",
                manifest=str(genotype_manifest),
                outputs={"genotypes_dir": str(genotypes_dir)},
            )
        )

        gene_lists_dir = out_dir / "gene_lists"
        gene_list_path, gene_list_manifest, gene_list_manifest_json = resolve_gene_source(args, gene_lists_dir)
        if gene_list_manifest is not None:
            steps.append(
                PipelineStepRecord(
                    name="build-gene-lists",
                    status="ok",
                    manifest=str(gene_list_manifest_json) if gene_list_manifest_json else None,
                    outputs={"gene_list_manifest": str(gene_list_manifest), "gene_lists_dir": str(gene_lists_dir)},
                )
            )
        else:
            steps.append(
                PipelineStepRecord(
                    name="build-gene-lists",
                    status="skipped",
                    message="Using caller-provided shared GeneList.txt",
                    outputs={"gene_list": str(gene_list_path)},
                )
            )

        build_dir = out_dir / "build"
        build_manifest = build_dir / "build_manifest.json"
        build_args = argparse.Namespace(
            annotated_dir=annotated_dir,
            genotypes_dir=genotypes_dir,
            gene_list=gene_list_path,
            gene_list_manifest=gene_list_manifest,
            out_dir=build_dir,
            samples=None,
            all_samples=True,
            chromosomes=args.chromosomes,
            dosage_threshold=args.dosage_threshold,
            qc_mode=args.qc_mode,
            force=args.force,
            dry_run=False,
            python=args.python,
            rscript=args.rscript,
            no_optimize=args.no_optimize,
            keep_temp_vcf=False,
            write_preview_csv=False,
            skip_chrom_gds=False,
            merge_per_sample=True,
            manifest_json=build_manifest,
            flatten_script=package_path("impact_snv.gds.flatten"),
            gds_writer=resource_script("favor_flat_to_seqarray_gds.R"),
        )
        rc = run_build_gds(build_args)
        if rc != 0:
            raise PipelineRunError("build-gds failed")
        steps.append(
            PipelineStepRecord(
                name="build-gds",
                status="ok",
                manifest=str(build_manifest),
                outputs={"build_dir": str(build_dir), "gds_merged": str(build_dir / "gds_merged")},
            )
        )

        qc_dir = out_dir / "qc"
        qc_args = argparse.Namespace(manifest=build_manifest, out_dir=qc_dir, qc_mode=args.qc_mode)
        rc = run_qc_build(qc_args)
        if rc != 0 and args.qc_mode == "strict":
            raise PipelineRunError("qc-build failed in strict mode")
        steps.append(
            PipelineStepRecord(
                name="qc-build",
                status="ok" if rc == 0 else "warning",
                outputs={"qc_dir": str(qc_dir), "build_qc_summary": str(qc_dir / "build_qc_summary.json")},
            )
        )

        final_dir = out_dir / "final"
        finalize_manifest = final_dir / "finalize_manifest.json"
        finalize_args = argparse.Namespace(
            input_dir=build_dir / "gds_merged",
            gene_list=gene_list_path,
            gene_list_manifest=gene_list_manifest,
            out_dir=final_dir,
            qc_mode=args.qc_mode,
            samples=None,
            rscript=args.rscript,
            manifest_json=finalize_manifest,
            force=args.force,
            no_optimize=args.no_optimize,
            keep_intermediate=False,
        )
        rc = run_finalize_gds(finalize_args)
        if rc != 0:
            raise PipelineRunError("finalize-gds failed")
        steps.append(
            PipelineStepRecord(
                name="finalize-gds",
                status="ok",
                manifest=str(finalize_manifest),
                outputs={"final_dir": str(final_dir)},
            )
        )

        validation_dir = out_dir / "validate"
        validation_dir.mkdir(parents=True, exist_ok=True)
        validation_manifest = validation_dir / "validate_manifest.json"
        validation_results = validate_gds_dir(final_dir, rscript=args.rscript, qc_mode=args.qc_mode)
        write_validation_manifest(validation_results, validation_manifest)
        if any(not result.is_valid for result in validation_results):
            raise PipelineRunError("validate-gds failed")
        steps.append(
            PipelineStepRecord(
                name="validate-gds",
                status="ok",
                manifest=str(validation_manifest),
                outputs={"validation_manifest": str(validation_manifest)},
            )
        )

        if not args.quiet:
            print("impact-snv run completed successfully")
            print(f"  cases: {len(case_ids)}")
            print(f"  merged_vcf: {merged_vcf}")
            print(f"  annotated_dir: {annotated_dir}")
            print(f"  final_dir: {final_dir}")
            print(f"  run_manifest: {manifest_path}")

    except Exception as exc:
        status = "failed"
        message = str(exc)
        print(f"ERROR: {message}")
    finally:
        payload = {
            "command": "impact-snv run",
            "version": VERSION,
            "created_utc": utc_now_iso(),
            "status": status,
            "message": message,
            "reference_build": reference_build,
            "samples_manifest": str(samples_manifest),
            "phenotypes_manifest": str(phenotypes_manifest) if phenotypes_manifest else None,
            "gene_list": str(Path(args.gene_list).expanduser().resolve()) if getattr(args, "gene_list", None) else None,
            "out_dir": str(out_dir),
            "run_prefix": args.run_prefix,
            "sanitize_inputs": args.sanitize_inputs,
            "steps": [asdict(step) for step in steps],
        }
        write_run_manifest(manifest_path, payload)

    return 0 if status == "ok" else 1


__all__ = ["run_pipeline", "detect_duplicate_input_sample_ids", "collect_input_vcfs"]