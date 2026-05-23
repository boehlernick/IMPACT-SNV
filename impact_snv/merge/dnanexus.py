#!/usr/bin/env python3
"""DNAnexus wrapper for IMPACT-SNV VCF merge.

This wrapper intentionally delegates all platform-independent merge behavior to
impact_snv.merge.core so local and DNAnexus behavior remain aligned.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

import dxpy

from impact_snv.merge.core import MergeConfig, VCFMergeError, merge_result_manifest, merge_vcf_samples


def log(msg: str) -> None:
    try:
        dxpy.log(msg)
    except Exception:
        pass
    print(msg, flush=True)


def file_id(x: Any) -> str:
    if isinstance(x, str):
        return x
    if isinstance(x, dict):
        link = x.get("$dnanexus_link")
        if isinstance(link, str):
            return link
        if isinstance(link, dict) and "id" in link:
            return str(link["id"])
        if "id" in x:
            return str(x["id"])
    raise dxpy.AppError(f"Could not interpret DNAnexus file input: {x!r}")


def safe_name(name: Optional[str], fallback: str) -> str:
    n = Path(name or fallback).name
    return n if n not in {"", ".", ".."} else fallback


def download(x: Any, dest: Path, fallback: str) -> Path:
    fid = file_id(x)
    desc = dxpy.DXFile(fid).describe()
    p = dest / safe_name(desc.get("name"), fallback)
    if p.exists():
        p = dest / f"{fid}_{p.name}"
    log(f"Downloading {fid} -> {p}")
    dxpy.download_dxfile(fid, str(p))
    return p


def upload(path: Optional[Path], name: Optional[str] = None):
    if not path:
        return None
    p = Path(path)
    if not p.exists() or not p.is_file():
        return None
    dx = dxpy.upload_local_file(str(p), name=name or p.name)
    return dxpy.dxlink(dx)


def out_name(name: Optional[str]) -> str:
    n = Path(name or "merged.vcf.gz").name
    if n.endswith(".vcf.gz"):
        return n
    if n.endswith(".vcf"):
        return n + ".gz"
    if n.endswith(".gz"):
        return n
    return n + ".vcf.gz"


def mode_from(normalization_mode: Optional[str], normalize: Optional[bool]) -> str:
    if normalization_mode is not None:
        mode = str(normalization_mode).lower()
    elif normalize is True:
        mode = "always"
    elif normalize is False:
        mode = "never"
    else:
        mode = "auto"
    if mode not in {"auto", "always", "never"}:
        raise dxpy.AppError("normalization_mode must be auto, always, or never")
    return mode


def upload_intermediates(paths: Iterable[Path], keep: bool) -> list[dict]:
    results = []
    if not keep:
        return results
    for p in paths:
        for q in [p, Path(str(p) + ".tbi"), Path(str(p) + ".csi")]:
            u = upload(q)
            if u:
                results.append(u)
    return results


@dxpy.entry_point("main")
def main(
    input_vcfs: Sequence[Any],
    reference_genome: Optional[Any] = None,
    reference_fai: Optional[Any] = None,
    output_name: str = "merged.vcf.gz",
    normalization_mode: Optional[str] = "auto",
    normalize: Optional[bool] = None,
    preflight_records: int = 10000,
    threads: int = 4,
    force_samples: bool = False,
    keep_intermediates: bool = False,
    extra_merge_args: Optional[Sequence[str]] = None,
    reference_build: str = "GRCh38",
    qc_mode: str = "warn",
    **kwargs: Any,
) -> dict[str, Any]:
    if reference_build != "GRCh38":
        raise dxpy.AppError("IMPACT-SNV v1.0.0 supports only reference_build=GRCh38")
    if not input_vcfs or len(input_vcfs) < 2:
        raise dxpy.AppError("input_vcfs must contain at least two VCF files")

    mode = mode_from(normalization_mode, normalize)
    if mode == "always" and reference_genome is None:
        raise dxpy.AppError("reference_genome is required when normalization_mode='always'")

    root = Path.cwd().resolve()
    inputs = root / "inputs"
    work = root / "work"
    outputs = root / "outputs"
    for d in (inputs, work, outputs):
        d.mkdir(parents=True, exist_ok=True)

    local_vcfs = [download(v, inputs, f"input_{i}.vcf.gz") for i, v in enumerate(input_vcfs, 1)]
    local_ref = None
    if reference_genome is not None and mode != "never":
        local_ref = download(reference_genome, inputs, "reference.fa")
        if reference_fai is not None:
            got = download(reference_fai, inputs, local_ref.name + ".fai")
            expected = Path(str(local_ref) + ".fai")
            if got != expected:
                if expected.exists():
                    expected.unlink()
                got.rename(expected)
            log(f"Reference FASTA index available at {expected}")

    config = MergeConfig(
        reference_fasta=local_ref,
        normalization_mode=mode,
        preflight_records=int(preflight_records or 10000),
        threads=int(threads or 1),
        keep_intermediates=bool(keep_intermediates),
        force_samples=bool(force_samples),
        overwrite=True,
        work_dir=work,
        extra_merge_args=tuple(str(x) for x in (extra_merge_args or [])),
    )

    try:
        result = merge_vcf_samples(local_vcfs, outputs / out_name(output_name), config)
    except VCFMergeError as e:
        raise dxpy.AppError(str(e)) from e

    cmd_log = outputs / "merge_command_log.txt"
    cmd_log.write_text("\n".join(result.command_log) + "\n", encoding="utf-8")
    manifest = {
        "command": "impact-snv merge",
        "status": "ok",
        **merge_result_manifest(result, qc_mode=qc_mode, reference_build=reference_build),
    }
    manifest_path = outputs / "merge_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    return {
        "merged_vcf": upload(result.merged_vcf, out_name(output_name)),
        "merged_vcf_index": upload(result.merged_index),
        "merge_command_log": upload(cmd_log),
        "merge_manifest_json": upload(manifest_path),
        "intermediate_vcfs": upload_intermediates(result.prepared_vcfs, bool(keep_intermediates)),
    }


dxpy.run()
