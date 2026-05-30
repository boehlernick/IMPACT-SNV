#!/usr/bin/env python3
"""FAVOR CLI annotate wrapper for IMPACT-SNV.

Annotation-only success is valid for current FAVOR CLI. Genotype parquet output
is produced by impact-snv extract-genotypes for v1.0.0.
"""
from __future__ import annotations
from pathlib import Path
from typing import Any
from impact_snv.favor.common import FavorRunResult, chromosome_partition_dirs, ensure_out_dir, observed_mtimes, observed_sizes, remove_existing_outputs, resolve_favor_bin, run_logged_command_with_progress, utc_now_iso, validate_reference_build, write_manifest

def infer_prefix_from_ingested_dir(path: Path) -> str:
    return path.name[:-len('.ingested')] if path.name.endswith('.ingested') else path.name

def _warning(code: str, message: str) -> dict[str, str]:
    return {"code": code, "severity": "warning", "message": message}

def run_favor_annotate(args: Any) -> int:
    ingested_dir=Path(args.ingested_dir).expanduser().resolve()
    if not ingested_dir.exists() or not ingested_dir.is_dir(): raise FileNotFoundError(f"Ingested directory not found: {ingested_dir}")
    reference_build=getattr(args,'reference_build','GRCh38'); validate_reference_build(reference_build)
    favor_bin=resolve_favor_bin(getattr(args,'favor_bin','favor')); out_dir=ensure_out_dir(Path(args.out_dir)); out_prefix=getattr(args,'out_prefix',None) or infer_prefix_from_ingested_dir(ingested_dir); force=bool(getattr(args,'force',False))
    progress_interval=int(getattr(args,'progress_interval_seconds',60) or 60); quiet=bool(getattr(args,'quiet',False)); tail_log_lines=int(getattr(args,'tail_log_lines',0) or 0); max_chars=int(getattr(args,'max_progress_log_line_chars',300) or 300); progress_mode=getattr(args,'progress_mode','normal') or 'normal'
    annotated_dir=out_dir/f"{out_prefix}.annotated"; genotypes_dir=out_dir/f"{out_prefix}.genotypes"; samples_file=genotypes_dir/'samples.txt'; removed_paths=remove_existing_outputs([annotated_dir], force=force)
    try: ingested_arg=str(ingested_dir.relative_to(out_dir))
    except ValueError: ingested_arg=str(ingested_dir)
    stdout_log=out_dir/'logs'/f"{out_prefix}.favor_annotate.stdout.log"; stderr_log=out_dir/'logs'/f"{out_prefix}.favor_annotate.stderr.log"; command=[favor_bin,'annotate',ingested_arg]
    start=utc_now_iso(); rc=run_logged_command_with_progress(command,cwd=out_dir,stdout_log=stdout_log,stderr_log=stderr_log,progress_label='favor-annotate',progress_interval_seconds=progress_interval,quiet=quiet,watch_paths=[annotated_dir],tail_log_lines=tail_log_lines,max_progress_log_line_chars=max_chars,progress_mode=progress_mode); end=utc_now_iso()
    ann_parts=chromosome_partition_dirs(annotated_dir); gt_parts=chromosome_partition_dirs(genotypes_dir); status='ok'; messages=[]; warnings=[]
    if rc!=0: status='failed'; messages.append(f"favor annotate exited with return code {rc}; see {stderr_log}")
    if not annotated_dir.exists(): status='failed'; messages.append(f"Expected annotated output was not created: {annotated_dir}")
    if not ann_parts: status='failed'; messages.append('No chromosome partition directories were found in annotated output')
    if not genotypes_dir.exists(): warnings.append(_warning('GENOTYPES_OUTPUT_NOT_CREATED', f"favor annotate did not create {genotypes_dir}; genotype extraction is handled by impact-snv extract-genotypes."))
    if not samples_file.exists(): warnings.append(_warning('GENOTYPE_SAMPLES_FILE_NOT_CREATED', f"favor annotate did not create {samples_file}; genotype extraction is handled by impact-snv extract-genotypes."))
    msg='; '.join(messages)
    result=FavorRunResult(command='favor-annotate',favor_command=list(map(str,command)),cwd=str(out_dir),out_dir=str(out_dir),out_prefix=out_prefix,reference_build=reference_build,favor_bin=getattr(args,'favor_bin','favor'),force=force,return_code=rc,start_time_utc=start,end_time_utc=end,stdout_log=str(stdout_log),stderr_log=str(stderr_log),progress_interval_seconds=progress_interval,quiet=quiet,tail_log_lines=tail_log_lines,max_progress_log_line_chars=max_chars,progress_mode=progress_mode,removed_paths=removed_paths,observed_output_sizes_bytes=observed_sizes({'annotated_dir':annotated_dir}),observed_output_latest_mtime_utc=observed_mtimes({'annotated_dir':annotated_dir}),status=status,message=msg)
    payload=result.to_dict(); payload.update({'ingested_dir':str(ingested_dir),'annotated_dir':str(annotated_dir.resolve()),'genotypes_dir':str(genotypes_dir.resolve()),'samples_file':str(samples_file.resolve()),'annotated_chromosome_partitions':ann_parts,'genotype_chromosome_partitions':gt_parts,'warnings':warnings})
    manifest=Path(getattr(args,'manifest_json',None) or (out_dir/f"{out_prefix}.annotate_manifest.json")); write_manifest(manifest,payload)
    if status!='ok': print(msg); return 1
    if warnings:
        for w in warnings: print(f"WARNING [{w['code']}]: {w['message']}")
    print(f"FAVOR annotate completed successfully: {annotated_dir}"); print(f"Manifest written: {manifest}"); return 0
