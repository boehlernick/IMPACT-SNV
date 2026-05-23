#!/usr/bin/env python3
"""Core VCF merge utilities for IMPACT step1_vcf_merge.

This module is platform independent: no DNAnexus and no argparse. It prepares
VCF/VCF.GZ/BCF inputs, performs heuristic auto-normalization preflight by
default, optionally runs bcftools norm, merges samples, and indexes the result.

Auto-normalization warning: the preflight is a fast heuristic. It is not the
same as full bcftools normalization against the original FASTA and cannot prove
that a VCF is normalized. It only detects obvious issues such as multiallelic
ALT alleles or non-parsimonious simple indels in inspected records.
"""
from __future__ import annotations

import os, shutil, subprocess, tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

VALID_NORMALIZATION_MODES = {"auto", "always", "never"}

class VCFMergeError(RuntimeError):
    """Raised when the merge pipeline cannot complete safely."""

@dataclass(frozen=True)
class NormalizationPreflight:
    mode: str
    should_normalize: bool
    inspected_records: int
    reasons: Tuple[str, ...]
    warning: str

@dataclass(frozen=True)
class MergeConfig:
    reference_fasta: Optional[Path] = None
    normalization_mode: str = "auto"  # auto | always | never
    preflight_records: int = 10000
    threads: int = 1
    keep_intermediates: bool = False
    force_samples: bool = False
    overwrite: bool = False
    work_dir: Optional[Path] = None
    bcftools_path: Optional[Path] = None
    bgzip_path: Optional[Path] = None
    extra_merge_args: Sequence[str] = field(default_factory=tuple)

@dataclass(frozen=True)
class MergeResult:
    merged_vcf: Path
    merged_index: Optional[Path]
    input_vcfs: Tuple[Path, ...]
    prepared_vcfs: Tuple[Path, ...]
    sample_ids: Tuple[str, ...]
    command_log: Tuple[str, ...]
    normalization_preflight: NormalizationPreflight

def resolve_executable(name: str, explicit_path: Optional[Path] = None) -> Path:
    if explicit_path:
        p = Path(explicit_path).expanduser().resolve()
        if p.exists() and os.access(p, os.X_OK): return p
        raise VCFMergeError(f"Executable for {name!r} is not executable: {p}")
    found = shutil.which(name)
    if found: return Path(found)
    raise VCFMergeError(f"Required executable {name!r} was not found on PATH")

def run_command(cmd: Sequence[os.PathLike|str], *, description: str, command_log: Optional[List[str]]=None, capture_stdout: bool=False):
    cmd_s = " ".join(str(x) for x in cmd)
    if command_log is not None: command_log.append(cmd_s)
    try:
        if capture_stdout:
            return subprocess.run([str(x) for x in cmd], check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return subprocess.run([str(x) for x in cmd], check=True)
    except subprocess.CalledProcessError as e:
        stderr = f"\nStderr:\n{e.stderr}" if getattr(e, 'stderr', None) else ""
        raise VCFMergeError(f"Command failed while attempting to {description}.\nExit code: {e.returncode}\nCommand: {cmd_s}{stderr}") from e

def assert_readable_file(path: Path, label="file") -> Path:
    p = Path(path).expanduser().resolve()
    if not p.exists(): raise VCFMergeError(f"Input {label} does not exist: {p}")
    if not p.is_file(): raise VCFMergeError(f"Input {label} is not a regular file: {p}")
    if not os.access(p, os.R_OK): raise VCFMergeError(f"Input {label} is not readable: {p}")
    return p

def assert_reference_ready(reference_fasta: Path) -> Path:
    ref = assert_readable_file(reference_fasta, "reference FASTA")
    fai = Path(str(ref)+".fai")
    if not fai.exists(): raise VCFMergeError(f"Reference FASTA index is missing: {fai}\nCreate it with: samtools faidx <reference.fa>")
    return ref

def validate_config(config: MergeConfig) -> MergeConfig:
    if config.normalization_mode not in VALID_NORMALIZATION_MODES:
        raise VCFMergeError(f"normalization_mode must be one of {sorted(VALID_NORMALIZATION_MODES)}")
    if config.threads < 1: raise VCFMergeError("threads must be >= 1")
    if config.preflight_records < 1: raise VCFMergeError("preflight_records must be >= 1")
    if config.normalization_mode == "always" and config.reference_fasta is None:
        raise VCFMergeError("normalization_mode='always' requires reference_fasta")
    if config.reference_fasta is not None:
        object.__setattr__(config, "reference_fasta", assert_reference_ready(config.reference_fasta))
    return config

def normalize_output_vcf_path(output_vcf: os.PathLike|str) -> Path:
    p = Path(output_vcf).expanduser().resolve(); n = p.name
    if n.endswith((".vcf.gz", ".bcf", ".gz")): return p
    if n.endswith(".vcf"): return p.with_name(n+".gz")
    return p.with_name(n+".vcf.gz")

def existing_vcf_index(vcf: Path) -> Optional[Path]:
    for s in (".tbi", ".csi"):
        idx = Path(str(vcf)+s)
        if idx.exists(): return idx
    return None

def ensure_bgzipped(vcf: Path, *, bgzip: Path, work_dir: Path, command_log: List[str], overwrite=False) -> Path:
    vcf = assert_readable_file(vcf, "VCF")
    if vcf.name.endswith((".vcf.gz", ".bcf")): return vcf
    if not vcf.name.endswith(".vcf"): raise VCFMergeError(f"Unsupported input extension: {vcf}")
    out = work_dir / f"{vcf.name}.gz"
    if out.exists() and not overwrite: raise VCFMergeError(f"Refusing to overwrite intermediate: {out}")
    cmd=[bgzip,"-c",vcf]; command_log.append(" ".join(str(x) for x in cmd)+f" > {out}")
    with out.open("wb") as fh: subprocess.run([str(x) for x in cmd], check=True, stdout=fh)
    return out

def index_vcf(vcf: Path, *, bcftools: Path, threads: int, command_log: List[str], force=False) -> Optional[Path]:
    vcf = assert_readable_file(vcf, "VCF to index")
    old = existing_vcf_index(vcf)
    if old and not force: return old
    cmd=[bcftools,"index","--threads",str(max(1,threads))]
    if force: cmd.append("--force")
    cmd.append(vcf)
    run_command(cmd, description=f"index {vcf}", command_log=command_log)
    return existing_vcf_index(vcf)

def validate_vcf(vcf: Path, *, bcftools: Path, command_log: List[str]) -> None:
    run_command([bcftools,"view","--no-version","--header-only",vcf], description=f"validate header {vcf}", command_log=command_log, capture_stdout=True)

def get_sample_ids(vcf: Path, *, bcftools: Path, command_log: List[str]) -> List[str]:
    r = run_command([bcftools,"query","-l",vcf], description=f"read sample IDs from {vcf}", command_log=command_log, capture_stdout=True)
    return [x.strip() for x in r.stdout.splitlines() if x.strip()]

def assert_unique_sample_ids(sample_ids_by_file: Mapping[Path, Sequence[str]]) -> Tuple[str, ...]:
    seen: Dict[str, Path] = {}; dup=[]; ordered=[]
    for p, ids in sample_ids_by_file.items():
        for sid in ids:
            ordered.append(sid)
            if sid in seen: dup.append((sid, seen[sid], p))
            else: seen[sid]=p
    if dup:
        msg="\n".join(f"  sample {s!r}: {a} and {b}" for s,a,b in dup)
        raise VCFMergeError("Duplicate sample IDs were detected across input VCFs. Reheader samples before merging.\nDuplicates:\n"+msg+"\nIf intentional, set force_samples=True.")
    return tuple(ordered)

def _record_issue(fields: Sequence[str]) -> Optional[str]:
    if len(fields)<5: return "malformed VCF record with fewer than 5 columns"
    ref, alt_field = fields[3], fields[4]
    if "," in alt_field: return "multiallelic ALT allele detected"
    if alt_field in {".","*"} or alt_field.startswith("<") or "]" in alt_field or "[" in alt_field: return None
    for alt in alt_field.split(','):
        if len(ref)>1 and len(alt)>1:
            if ref[-1].upper()==alt[-1].upper(): return "non-parsimonious REF/ALT sharing a suffix"
            if len(ref)>2 and len(alt)>2 and ref[0].upper()==alt[0].upper(): return "potentially non-parsimonious REF/ALT sharing a long prefix"
    return None

def preflight_normalization_check(vcf_paths: Sequence[Path], *, bcftools: Path, max_records_per_file: int, command_log: List[str]) -> NormalizationPreflight:
    warning = "Auto-normalization uses a lightweight heuristic preflight. It is not equivalent to full bcftools normalization against the original reference FASTA and cannot prove that a VCF is fully normalized."
    reasons=[]; inspected=0
    for vcf in vcf_paths:
        cmd=[str(bcftools),"view","-H",str(vcf)]; command_log.append(" ".join(cmd)+f" | head -n {max_records_per_file}")
        proc=subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        assert proc.stdout is not None
        stderr=""
        try:
            for i,line in enumerate(proc.stdout,1):
                if i>max_records_per_file: break
                inspected += 1
                issue=_record_issue(line.rstrip('\n').split('\t'))
                if issue:
                    reasons.append(f"{vcf}: {issue} near inspected record {i}"); break
        finally:
            proc.stdout.close();
            try: proc.kill()
            except OSError: pass
            stderr = proc.stderr.read() if proc.stderr else ""; proc.wait()
            if proc.stderr: proc.stderr.close()
        if proc.returncode not in (0,-9,-13,141) and stderr.strip(): raise VCFMergeError(f"bcftools preflight failed for {vcf}:\n{stderr}")
    return NormalizationPreflight("auto", bool(reasons), inspected, tuple(reasons), warning)

def normalize_and_split_vcf(input_vcf: Path, *, output_vcf: Path, reference_fasta: Path, bcftools: Path, threads: int, command_log: List[str], overwrite=False) -> Path:
    if output_vcf.exists() and not overwrite: raise VCFMergeError(f"Refusing to overwrite normalized VCF: {output_vcf}")
    run_command([bcftools,"norm","--threads",str(max(1,threads)),"-m","-any","-f",reference_fasta,"-Oz","-o",output_vcf,input_vcf], description=f"normalize {input_vcf}", command_log=command_log)
    return output_vcf

def decide_normalization(prepared: Sequence[Path], *, config: MergeConfig, bcftools: Path, command_log: List[str]) -> NormalizationPreflight:
    if config.normalization_mode == "never": return NormalizationPreflight("never", False, 0, ("Normalization explicitly disabled by user",), "Normalization was skipped by request; VCF records are merged as-is.")
    if config.normalization_mode == "always": return NormalizationPreflight("always", True, 0, ("Normalization explicitly requested by user",), "Full bcftools norm will be run for all prepared input VCFs.")
    pf=preflight_normalization_check(prepared, bcftools=bcftools, max_records_per_file=config.preflight_records, command_log=command_log)
    if pf.should_normalize and config.reference_fasta is None:
        raise VCFMergeError("Auto-normalization preflight found records that appear to need normalization, but no reference FASTA was provided. Provide --reference or use --no-normalize.\n"+"\n".join(f"  - {r}" for r in pf.reasons)+"\n\n"+pf.warning)
    return pf

def merge_vcf_samples(input_vcfs: Sequence[os.PathLike|str], output_vcf: os.PathLike|str, config: Optional[MergeConfig]=None) -> MergeResult:
    config = validate_config(config or MergeConfig()); log: List[str]=[]
    if len(input_vcfs)<2: raise VCFMergeError("At least two input VCFs are required")
    inputs=tuple(assert_readable_file(Path(p), "input VCF") for p in input_vcfs)
    bcftools=resolve_executable("bcftools", config.bcftools_path); bgzip=resolve_executable("bgzip", config.bgzip_path)
    tmp=None
    if config.work_dir: work=Path(config.work_dir).expanduser().resolve(); work.mkdir(parents=True, exist_ok=True)
    elif config.keep_intermediates: work=Path(tempfile.mkdtemp(prefix="vcf_merge_"))
    else: tmp=tempfile.TemporaryDirectory(prefix="vcf_merge_"); work=Path(tmp.name)
    try:
        basic=[]
        for inp in inputs:
            v=ensure_bgzipped(inp, bgzip=bgzip, work_dir=work, command_log=log, overwrite=config.overwrite)
            validate_vcf(v, bcftools=bcftools, command_log=log); index_vcf(v, bcftools=bcftools, threads=config.threads, command_log=log)
            basic.append(v)
        pf=decide_normalization(basic, config=config, bcftools=bcftools, command_log=log)
        prepared=[]
        if pf.should_normalize:
            for v in basic:
                stem=v.name[:-7] if v.name.endswith('.vcf.gz') else (v.name[:-4] if v.name.endswith('.bcf') else v.name)
                nv=normalize_and_split_vcf(v, output_vcf=work/f"{stem}.norm.vcf.gz", reference_fasta=config.reference_fasta, bcftools=bcftools, threads=config.threads, command_log=log, overwrite=config.overwrite)
                index_vcf(nv, bcftools=bcftools, threads=config.threads, command_log=log); prepared.append(nv)
        else: prepared=list(basic)
        sid_by_file={v:get_sample_ids(v, bcftools=bcftools, command_log=log) for v in prepared}
        sample_ids=tuple(s for ids in sid_by_file.values() for s in ids) if config.force_samples else assert_unique_sample_ids(sid_by_file)
        out=normalize_output_vcf_path(output_vcf); out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists() and not config.overwrite: raise VCFMergeError(f"Refusing to overwrite existing output VCF: {out}")
        cmd=[bcftools,"merge","--threads",str(max(1,config.threads)),"-Oz","-o",out]
        if config.force_samples: cmd.append("--force-samples")
        cmd.extend(config.extra_merge_args); cmd.extend(prepared)
        run_command(cmd, description="merge prepared VCFs", command_log=log)
        idx=index_vcf(out, bcftools=bcftools, threads=config.threads, command_log=log, force=True)
        return MergeResult(out, idx, inputs, tuple(prepared), sample_ids, tuple(log), pf)
    finally:
        if tmp is not None and not config.keep_intermediates: tmp.cleanup()

__all__=["MergeConfig","MergeResult","NormalizationPreflight","VCFMergeError","merge_vcf_samples","normalize_output_vcf_path","preflight_normalization_check","index_vcf","get_sample_ids"]
