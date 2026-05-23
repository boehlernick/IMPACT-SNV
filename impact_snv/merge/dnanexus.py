#!/usr/bin/env python3
"""DNAnexus wrapper for IMPACT step1_vcf_merge."""
from __future__ import annotations
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
import dxpy
from merge_core import MergeConfig, VCFMergeError, merge_vcf_samples

def log(msg):
    try: dxpy.log(msg)
    except Exception: pass
    print(msg, flush=True)

def file_id(x: Any) -> str:
    if isinstance(x,str): return x
    if isinstance(x,dict):
        link=x.get('$dnanexus_link')
        if isinstance(link,str): return link
        if isinstance(link,dict) and 'id' in link: return str(link['id'])
        if 'id' in x: return str(x['id'])
    raise dxpy.AppError(f'Could not interpret DNAnexus file input: {x!r}')

def safe_name(name, fallback):
    n=Path(name or fallback).name
    return n if n not in {'','.','..'} else fallback

def download(x, dest: Path, fallback: str) -> Path:
    fid=file_id(x); desc=dxpy.DXFile(fid).describe(); p=dest/safe_name(desc.get('name'), fallback)
    if p.exists(): p=dest/f'{fid}_{p.name}'
    log(f'Downloading {fid} -> {p}'); dxpy.download_dxfile(fid, str(p)); return p

def upload(path: Optional[Path], name: Optional[str]=None):
    if not path: return None
    p=Path(path)
    if not p.exists() or not p.is_file(): return None
    dx=dxpy.upload_local_file(str(p), name=name or p.name); return dxpy.dxlink(dx)

def out_name(name: Optional[str]) -> str:
    n=Path(name or 'merged.vcf.gz').name
    if n.endswith('.vcf.gz'): return n
    if n.endswith('.vcf'): return n+'.gz'
    if n.endswith('.gz'): return n
    return n+'.vcf.gz'

def mode_from(normalization_mode: Optional[str], normalize: Optional[bool]) -> str:
    m = str(normalization_mode).lower() if normalization_mode is not None else ('always' if normalize else 'never' if normalize is not None else 'auto')
    if m not in {'auto','always','never'}: raise dxpy.AppError('normalization_mode must be auto, always, or never')
    return m

def upload_intermediates(paths: Iterable[Path], keep: bool) -> List[dict]:
    res=[]
    if not keep: return res
    for p in paths:
        for q in [p, Path(str(p)+'.tbi'), Path(str(p)+'.csi')]:
            u=upload(q)
            if u: res.append(u)
    return res

@dxpy.entry_point('main')
def main(input_vcfs: Sequence[Any], reference_genome: Optional[Any]=None, reference_fai: Optional[Any]=None, output_name: str='merged.vcf.gz', normalization_mode: Optional[str]='auto', normalize: Optional[bool]=None, preflight_records: int=10000, threads: int=4, force_samples: bool=False, keep_intermediates: bool=False, extra_merge_args: Optional[Sequence[str]]=None, **kwargs) -> Dict[str, Any]:
    if not input_vcfs or len(input_vcfs)<2: raise dxpy.AppError('input_vcfs must contain at least two VCF files')
    mode=mode_from(normalization_mode, normalize)
    if mode=='always' and reference_genome is None: raise dxpy.AppError("reference_genome is required when normalization_mode='always'")
    root=Path.cwd().resolve(); inputs=root/'inputs'; work=root/'work'; outputs=root/'outputs'
    for d in (inputs,work,outputs): d.mkdir(parents=True, exist_ok=True)
    log(f'normalization_mode={mode}; threads={threads}; force_samples={force_samples}')
    local_vcfs=[download(v, inputs, f'input_{i}.vcf.gz') for i,v in enumerate(input_vcfs,1)]
    local_ref=None
    if reference_genome is not None and mode!='never':
        local_ref=download(reference_genome, inputs, 'reference.fa')
        if reference_fai is not None:
            got=download(reference_fai, inputs, local_ref.name+'.fai'); expected=Path(str(local_ref)+'.fai')
            if got != expected:
                if expected.exists(): expected.unlink()
                got.rename(expected)
            log(f'Reference FASTA index available at {expected}')
    config=MergeConfig(reference_fasta=local_ref, normalization_mode=mode, preflight_records=int(preflight_records or 10000), threads=int(threads or 1), keep_intermediates=bool(keep_intermediates), force_samples=bool(force_samples), overwrite=True, work_dir=work, extra_merge_args=tuple(str(x) for x in (extra_merge_args or [])))
    try: result=merge_vcf_samples(local_vcfs, outputs/out_name(output_name), config)
    except VCFMergeError as e: raise dxpy.AppError(str(e)) from e
    pf=result.normalization_preflight; log(pf.warning); log(f'Normalization performed: {pf.should_normalize}')
    cmd_log=outputs/'merge_command_log.txt'; cmd_log.write_text('\n'.join(result.command_log)+'\n')
    summ={"merged_vcf":str(result.merged_vcf),"merged_index":str(result.merged_index) if result.merged_index else None,"input_vcfs":[str(x) for x in result.input_vcfs],"prepared_vcfs":[str(x) for x in result.prepared_vcfs],"sample_ids":list(result.sample_ids),"sample_count":len(result.sample_ids),"normalization_mode":pf.mode,"normalization_performed":pf.should_normalize,"normalization_preflight_records":pf.inspected_records,"normalization_reasons":list(pf.reasons),"normalization_warning":pf.warning,"threads":int(threads or 1),"force_samples":bool(force_samples),"extra_merge_args":list(extra_merge_args or [])}
    summ_path=outputs/'merge_summary.json'; summ_path.write_text(json.dumps(summ, indent=2, sort_keys=True))
    return {"merged_vcf":upload(result.merged_vcf, out_name(output_name)), "merged_vcf_index":upload(result.merged_index), "merge_command_log":upload(cmd_log), "merge_summary_json":upload(summ_path), "intermediate_vcfs":upload_intermediates(result.prepared_vcfs, bool(keep_intermediates))}

dxpy.run()
