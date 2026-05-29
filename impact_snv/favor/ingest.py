#!/usr/bin/env python3
"""FAVOR CLI ingest wrapper for IMPACT-SNV."""
from __future__ import annotations
from pathlib import Path
from typing import Any
from impact_snv.favor.common import FavorRunResult, ensure_out_dir, existing_vcf_index, infer_prefix_from_vcf, observed_mtimes, observed_sizes, remove_existing_outputs, resolve_favor_bin, run_logged_command_with_progress, utc_now_iso, validate_reference_build, write_manifest

def run_favor_ingest(args: Any) -> int:
    input_vcf=Path(args.input_vcf).expanduser().resolve()
    if not input_vcf.exists() or not input_vcf.is_file(): raise FileNotFoundError(f"Input VCF not found: {input_vcf}")
    input_index=existing_vcf_index(input_vcf)
    if input_index is None: raise FileNotFoundError(f"Input VCF index not found; expected {input_vcf}.csi or {input_vcf}.tbi")
    reference_build=getattr(args,'reference_build','GRCh38'); validate_reference_build(reference_build)
    favor_bin=resolve_favor_bin(getattr(args,'favor_bin','favor')); out_dir=ensure_out_dir(Path(args.out_dir)); out_prefix=getattr(args,'out_prefix',None) or infer_prefix_from_vcf(input_vcf); force=bool(getattr(args,'force',False))
    progress_interval=int(getattr(args,'progress_interval_seconds',60) or 60); quiet=bool(getattr(args,'quiet',False)); tail_log_lines=int(getattr(args,'tail_log_lines',0) or 0); max_chars=int(getattr(args,'max_progress_log_line_chars',300) or 300); progress_mode=getattr(args,'progress_mode','normal') or 'normal'
    ingested_dir=out_dir/f"{out_prefix}.ingested"; removed_paths=remove_existing_outputs([ingested_dir], force=force)
    stdout_log=out_dir/'logs'/f"{out_prefix}.favor_ingest.stdout.log"; stderr_log=out_dir/'logs'/f"{out_prefix}.favor_ingest.stderr.log"; command=[favor_bin,'ingest',str(input_vcf)]
    start=utc_now_iso(); rc=run_logged_command_with_progress(command,cwd=out_dir,stdout_log=stdout_log,stderr_log=stderr_log,progress_label='favor-ingest',progress_interval_seconds=progress_interval,quiet=quiet,watch_paths=[ingested_dir],tail_log_lines=tail_log_lines,max_progress_log_line_chars=max_chars,progress_mode=progress_mode); end=utc_now_iso()
    status='ok' if rc==0 and ingested_dir.exists() else 'failed'; msg='' if status=='ok' else (f"favor ingest exited with return code {rc}; see {stderr_log}" if rc!=0 else f"Expected ingested output was not created: {ingested_dir}")
    result=FavorRunResult(command='favor-ingest',favor_command=list(map(str,command)),cwd=str(out_dir),out_dir=str(out_dir),out_prefix=out_prefix,reference_build=reference_build,favor_bin=getattr(args,'favor_bin','favor'),force=force,return_code=rc,start_time_utc=start,end_time_utc=end,stdout_log=str(stdout_log),stderr_log=str(stderr_log),progress_interval_seconds=progress_interval,quiet=quiet,tail_log_lines=tail_log_lines,max_progress_log_line_chars=max_chars,progress_mode=progress_mode,removed_paths=removed_paths,observed_output_sizes_bytes=observed_sizes({'ingested_dir':ingested_dir}),observed_output_latest_mtime_utc=observed_mtimes({'ingested_dir':ingested_dir}),status=status,message=msg)
    payload=result.to_dict(); payload.update({'input_vcf':str(input_vcf),'input_index':str(input_index),'ingested_dir':str(ingested_dir.resolve())})
    manifest=Path(getattr(args,'manifest_json',None) or (out_dir/f"{out_prefix}.ingest_manifest.json")); write_manifest(manifest,payload)
    if status!='ok': print(msg); return 1
    print(f"FAVOR ingest completed successfully: {ingested_dir}"); print(f"Manifest written: {manifest}"); return 0
