#!/usr/bin/env python3
"""VCF input sanitization helpers for stable unique sample IDs.

This module backs the production CLI command:

    impact-snv sanitize-vcfs

It reheaders one or more input VCF/VCF.GZ/BCF files so that sample IDs are
stable and unique across the set before merge/build/finalize workflows.
"""

from __future__ import annotations

import json
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
import re
from typing import Any, Optional, Sequence

from impact_snv import __version__ as VERSION
from impact_snv.merge.core import (
    VCFMergeError,
    assert_readable_file,
    assert_unique_sample_ids,
    ensure_bgzipped,
    get_sample_ids,
    index_vcf,
    resolve_executable,
    run_command,
    validate_vcf,
)


DEFAULT_PREFIX_REGEX = r"(?i)(case\d+)"


@dataclass(frozen=True)
class SanitizedVCFRecord:
    input_vcf: str
    prepared_input_vcf: str
    output_vcf: str
    output_index: Optional[str]
    sample_prefix: str
    original_sample_ids: tuple[str, ...]
    sanitized_sample_ids: tuple[str, ...]


@dataclass(frozen=True)
class SanitizeVCFsResult:
    command: str
    version: str
    created_utc: str
    out_dir: str
    manifest_tsv: str
    manifest_json: str
    sample_count: int
    records: tuple[SanitizedVCFRecord, ...]
    command_log: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            **asdict(self),
            "records": [asdict(record) for record in self.records],
        }


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def strip_vcf_suffix(name: str) -> str:
    for suffix in (".vcf.gz", ".vcf", ".bcf", ".gz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def normalize_token(text: str, fallback: str = "sample") -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text).strip())
    token = re.sub(r"_+", "_", token).strip("._-")
    return token or fallback


def derive_sample_prefix(path: Path, prefix_regex: str = DEFAULT_PREFIX_REGEX) -> str:
    stem = strip_vcf_suffix(path.name)
    stem = re.sub(r"\.reheadered$", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"\.sanitized$", "", stem, flags=re.IGNORECASE)
    if prefix_regex:
        match = re.search(prefix_regex, stem)
        if match:
            if match.lastindex:
                return normalize_token(match.group(1), fallback="sample")
            return normalize_token(match.group(0), fallback="sample")
    return normalize_token(stem, fallback="sample")


def assign_unique_prefixes(vcfs: Sequence[Path], prefix_regex: str = DEFAULT_PREFIX_REGEX) -> dict[Path, str]:
    assigned: dict[Path, str] = {}
    used: set[str] = set()
    for vcf in vcfs:
        base = derive_sample_prefix(vcf, prefix_regex=prefix_regex)
        candidate = base
        suffix = 2
        while candidate in used:
            candidate = f"{base}_{suffix}"
            suffix += 1
        used.add(candidate)
        assigned[vcf] = candidate
    return assigned


def build_sanitized_sample_ids(sample_ids: Sequence[str], prefix: str, seen_global: set[str]) -> tuple[str, ...]:
    local_seen: set[str] = set()
    out: list[str] = []
    for original in sample_ids:
        normalized_sample = normalize_token(original, fallback="sample")
        if normalized_sample == prefix or normalized_sample.startswith(f"{prefix}_"):
            base = normalized_sample
        else:
            base = normalize_token(f"{prefix}_{normalized_sample}", fallback=prefix)
        candidate = base
        suffix = 2
        while candidate in seen_global or candidate in local_seen:
            candidate = f"{base}_{suffix}"
            suffix += 1
        out.append(candidate)
        local_seen.add(candidate)
        seen_global.add(candidate)
    return tuple(out)


def write_sample_map(path: Path, sample_ids: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(sample_ids) + "\n", encoding="utf-8")


def output_vcf_path(out_dir: Path, prefix: str) -> Path:
    return out_dir / f"{normalize_token(prefix, fallback='sample')}.sanitized.vcf.gz"


def write_manifest_tsv(path: Path, records: Sequence[SanitizedVCFRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "vcf_path\toriginal_vcf\tsample_prefix\toriginal_sample_ids\tsanitized_sample_ids",
    ]
    for record in records:
        lines.append(
            "\t".join(
                [
                    record.output_vcf,
                    record.input_vcf,
                    record.sample_prefix,
                    ";".join(record.original_sample_ids),
                    ";".join(record.sanitized_sample_ids),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest_json(path: Path, result: SanitizeVCFsResult) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result.to_dict(), indent=2, sort_keys=True), encoding="utf-8")


def run_sanitize_vcfs(args: Any) -> int:
    vcfs = [Path(x).expanduser().resolve() for x in getattr(args, "vcfs", [])]
    if not vcfs:
        raise VCFMergeError("sanitize-vcfs requires at least one input VCF")

    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_tsv = Path(getattr(args, "manifest_tsv", None) or (out_dir / "sanitized_vcfs.tsv")).expanduser().resolve()
    manifest_json = Path(getattr(args, "manifest_json", None) or (out_dir / "sanitize_manifest.json")).expanduser().resolve()
    command_log_path = Path(getattr(args, "command_log", None)).expanduser().resolve() if getattr(args, "command_log", None) else None
    force = bool(getattr(args, "force", False))

    bcftools = resolve_executable("bcftools", Path(args.bcftools) if getattr(args, "bcftools", None) else None)
    needs_bgzip = any(not str(vcf).endswith((".vcf.gz", ".bcf")) for vcf in vcfs)
    bgzip = resolve_executable("bgzip", Path(args.bgzip) if getattr(args, "bgzip", None) else None) if needs_bgzip else None

    command_log: list[str] = []
    prefix_regex = getattr(args, "prefix_regex", DEFAULT_PREFIX_REGEX)
    threads = int(getattr(args, "threads", 1) or 1)

    with tempfile.TemporaryDirectory(prefix="vcf_sanitize_") as tmp_name:
        work_dir = Path(tmp_name)
        prepared_inputs: list[Path] = []
        sample_ids_by_file: dict[Path, tuple[str, ...]] = {}
        original_inputs_by_prepared: dict[Path, Path] = {}

        for input_vcf in vcfs:
            readable = assert_readable_file(input_vcf, "VCF")
            prepared = readable
            if bgzip is not None:
                prepared = ensure_bgzipped(readable, bgzip=bgzip, work_dir=work_dir, command_log=command_log, overwrite=force)
            validate_vcf(prepared, bcftools=bcftools, command_log=command_log)
            sample_ids = tuple(get_sample_ids(prepared, bcftools=bcftools, command_log=command_log))
            if not sample_ids:
                raise VCFMergeError(f"No sample IDs found in {prepared}")
            prepared_inputs.append(prepared)
            sample_ids_by_file[prepared] = sample_ids
            original_inputs_by_prepared[prepared] = readable

        prefixes = assign_unique_prefixes(prepared_inputs, prefix_regex=prefix_regex)
        seen_global: set[str] = set()
        records: list[SanitizedVCFRecord] = []
        sanitized_ids_by_file: dict[Path, tuple[str, ...]] = {}

        for prepared in prepared_inputs:
            prefix = prefixes[prepared]
            new_sample_ids = build_sanitized_sample_ids(sample_ids_by_file[prepared], prefix, seen_global)
            sanitized_ids_by_file[prepared] = new_sample_ids

        assert_unique_sample_ids(sanitized_ids_by_file)

        for prepared in prepared_inputs:
            prefix = prefixes[prepared]
            output_vcf = output_vcf_path(out_dir, prefix)
            if output_vcf.exists() and not force:
                raise VCFMergeError(f"Refusing to overwrite existing sanitized VCF: {output_vcf}")
            sample_map = work_dir / f"{normalize_token(prefix)}.samples.txt"
            write_sample_map(sample_map, sanitized_ids_by_file[prepared])
            run_command(
                [bcftools, "reheader", "-s", sample_map, "-o", output_vcf, prepared],
                description=f"reheader sample IDs for {prepared}",
                command_log=command_log,
            )
            output_index = index_vcf(output_vcf, bcftools=bcftools, threads=threads, command_log=command_log, force=True)
            written_ids = tuple(get_sample_ids(output_vcf, bcftools=bcftools, command_log=command_log))
            if written_ids != sanitized_ids_by_file[prepared]:
                raise VCFMergeError(
                    f"Sanitized sample IDs did not round-trip for {output_vcf}. "
                    f"Expected {sanitized_ids_by_file[prepared]}, observed {written_ids}"
                )
            records.append(
                SanitizedVCFRecord(
                    input_vcf=str(original_inputs_by_prepared[prepared]),
                    prepared_input_vcf=str(prepared),
                    output_vcf=str(output_vcf),
                    output_index=str(output_index) if output_index else None,
                    sample_prefix=prefix,
                    original_sample_ids=sample_ids_by_file[prepared],
                    sanitized_sample_ids=written_ids,
                )
            )

    result = SanitizeVCFsResult(
        command="impact-snv sanitize-vcfs",
        version=VERSION,
        created_utc=utc_now_iso(),
        out_dir=str(out_dir),
        manifest_tsv=str(manifest_tsv),
        manifest_json=str(manifest_json),
        sample_count=sum(len(record.sanitized_sample_ids) for record in records),
        records=tuple(records),
        command_log=tuple(command_log),
    )
    write_manifest_tsv(manifest_tsv, records)
    write_manifest_json(manifest_json, result)
    if command_log_path is not None:
        command_log_path.parent.mkdir(parents=True, exist_ok=True)
        command_log_path.write_text("\n".join(command_log) + "\n", encoding="utf-8")

    if not bool(getattr(args, "quiet", False)):
        print("VCF sample sanitization completed successfully")
        print(f"  files: {len(records)}")
        print(f"  samples: {result.sample_count}")
        print(f"  out_dir: {out_dir}")
        print(f"  manifest_tsv: {manifest_tsv}")
        print(f"  manifest_json: {manifest_json}")
        for record in records:
            print(f"  {Path(record.output_vcf).name}: {', '.join(record.sanitized_sample_ids)}")

    return 0


__all__ = [
    "DEFAULT_PREFIX_REGEX",
    "SanitizedVCFRecord",
    "SanitizeVCFsResult",
    "assign_unique_prefixes",
    "build_sanitized_sample_ids",
    "derive_sample_prefix",
    "run_sanitize_vcfs",
]