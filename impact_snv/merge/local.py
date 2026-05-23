#!/usr/bin/env python3
"""Local CLI wrapper for IMPACT step1_vcf_merge."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
from typing import Any, Dict, Optional, Sequence
from merge_core import MergeConfig, MergeResult, VCFMergeError, merge_vcf_samples

def positive_int(v: str) -> int:
    try: x=int(v)
    except ValueError as e: raise argparse.ArgumentTypeError(f"Expected integer, got {v!r}") from e
    if x<1: raise argparse.ArgumentTypeError("Value must be >= 1")
    return x

def build_parser():
    p=argparse.ArgumentParser(description="Merge VCF/VCF.GZ/BCF files into an indexed multi-sample VCF.GZ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--vcfs', nargs='+', required=True, help='Input VCF/VCF.GZ/BCF files; at least two')
    p.add_argument('--output', required=True, help='Output merged VCF. .vcf.gz is recommended; .vcf will be written as .vcf.gz')
    p.add_argument('--reference', help='Reference FASTA for bcftools norm; requires matching .fai')
    g=p.add_mutually_exclusive_group()
    g.add_argument('--auto-normalize', dest='normalization_mode', action='store_const', const='auto', default='auto', help='Heuristic preflight, normalize only if obvious issues are found')
    g.add_argument('--normalize', dest='normalization_mode', action='store_const', const='always', help='Always run bcftools norm')
    g.add_argument('--no-normalize', dest='normalization_mode', action='store_const', const='never', help='Merge as-is')
    p.add_argument('--preflight-records', type=positive_int, default=10000)
    p.add_argument('--threads', type=positive_int, default=1)
    p.add_argument('--work-dir')
    p.add_argument('--keep-intermediates', action='store_true')
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('--force-samples', action='store_true', help='Allow duplicate sample IDs via bcftools --force-samples')
    p.add_argument('--bcftools')
    p.add_argument('--bgzip')
    p.add_argument('--extra-merge-arg', action='append', default=[])
    p.add_argument('--command-log')
    p.add_argument('--summary-json')
    p.add_argument('--quiet', action='store_true')
    return p

def validate(args, parser):
    if len(args.vcfs)<2: parser.error('--vcfs requires at least two input files')
    if args.normalization_mode=='always' and not args.reference: parser.error('--normalize requires --reference')
    if args.normalization_mode=='never' and args.reference: parser.error('--reference is not used with --no-normalize')

def make_config(args):
    return MergeConfig(reference_fasta=Path(args.reference) if args.reference else None, normalization_mode=args.normalization_mode, preflight_records=args.preflight_records, threads=args.threads, keep_intermediates=args.keep_intermediates, force_samples=args.force_samples, overwrite=args.overwrite, work_dir=Path(args.work_dir) if args.work_dir else None, bcftools_path=Path(args.bcftools) if args.bcftools else None, bgzip_path=Path(args.bgzip) if args.bgzip else None, extra_merge_args=tuple(args.extra_merge_arg or []))

def summary(result: MergeResult, args) -> Dict[str, Any]:
    pf=result.normalization_preflight
    return {"merged_vcf":str(result.merged_vcf),"merged_index":str(result.merged_index) if result.merged_index else None,"input_vcfs":[str(x) for x in result.input_vcfs],"prepared_vcfs":[str(x) for x in result.prepared_vcfs],"sample_ids":list(result.sample_ids),"sample_count":len(result.sample_ids),"normalization_mode":pf.mode,"normalization_performed":pf.should_normalize,"normalization_preflight_records":pf.inspected_records,"normalization_reasons":list(pf.reasons),"normalization_warning":pf.warning,"threads":args.threads,"extra_merge_args":list(args.extra_merge_arg or []),"command_count":len(result.command_log)}

def write(path, text):
    p=Path(path).expanduser().resolve(); p.parent.mkdir(parents=True, exist_ok=True); p.write_text(text, encoding='utf-8')

def main(argv: Optional[Sequence[str]]=None) -> int:
    parser=build_parser(); args=parser.parse_args(argv); validate(args, parser)
    try: result=merge_vcf_samples([Path(x) for x in args.vcfs], Path(args.output), make_config(args))
    except VCFMergeError as e: print(f"ERROR: {e}", file=sys.stderr); return 1
    except KeyboardInterrupt: print('ERROR: interrupted by user', file=sys.stderr); return 130
    s=summary(result,args)
    if args.command_log: write(args.command_log, '\n'.join(result.command_log)+'\n')
    if args.summary_json: write(args.summary_json, json.dumps(s, indent=2, sort_keys=True))
    if not args.quiet:
        print('VCF merge completed successfully')
        print(f"  merged_vcf:              {result.merged_vcf}")
        print(f"  merged_index:            {result.merged_index or 'not found'}")
        print(f"  samples:                 {len(result.sample_ids)}")
        print(f"  normalization_mode:      {s['normalization_mode']}")
        print(f"  normalization_performed: {s['normalization_performed']}")
        print(f"  preflight_records:       {s['normalization_preflight_records']}")
        print(f"  warning:                 {s['normalization_warning']}")
        for r in s['normalization_reasons']: print(f"    - {r}")
    return 0
if __name__=='__main__': raise SystemExit(main())
