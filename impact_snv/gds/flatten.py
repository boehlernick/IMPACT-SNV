#!/usr/bin/env python3
"""
flatten_favor_for_gds.py

Flatten FAVOR annotated parquet + FAVOR genotype parquet into a compact,
per-sample, per-chromosome table suitable for validation before writing
SeqArray-compatible IMPACT GDS files.

This script implements the first-pass IMPACT GDS inclusion rule:

    include variant if:
      1. the selected sample has dosage > threshold, and
      2. the FAVOR annotation maps the variant to at least one gene in GeneList.txt

The output is NOT a GDS file. It is a flat preview/intermediate table that lets us
validate joins, gene filtering, dosage extraction, and annotation mapping before
we implement the R GDS writer.

Expected inputs
---------------
annotated_dir/
  chromosome=1/data.parquet
  chromosome=2/data.parquet
  ...

genotypes_dir/
  chromosome=1/data.parquet
  chromosome=2/data.parquet
  ...
  samples.txt

GeneList.txt:
  symbol<TAB>globalScore

Example
-------
python flatten_favor_for_gds.py \
  --annotated-dir /path/merged_case1_case2_output.annotated \
  --genotypes-dir /path/merged_case1_case2_output.genotypes \
  --gene-list /path/GeneList.txt \
  --sample-id Case1_proband \
  --chromosome 1 \
  --out Case1_proband_chr1.preview.parquet \
  --preview-csv Case1_proband_chr1.preview.head.csv \
  --summary-json Case1_proband_chr1.preview.summary.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
import pyarrow.parquet as pq


# -----------------------------
# Basic helpers
# -----------------------------

def normalize_chrom(value: Any) -> str:
    """Normalize chromosome labels to 1..22, X, Y, M/MT without chr prefix."""
    if value is None:
        return ""
    text = str(value)
    if text.lower().startswith("chr"):
        text = text[3:]
    # Avoid rendering numeric chromosomes as "1.0"
    if text.endswith(".0"):
        text = text[:-2]
    return text


def scalar_missing(value: Any) -> bool:
    """Robust missing-value test for scalar-ish values."""
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except Exception:
        return False


def as_list(value: Any) -> List[Any]:
    """Convert Arrow/Pandas list-like values to a normal Python list."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    # numpy arrays, pandas arrays, pyarrow list-like objects
    if hasattr(value, "tolist"):
        out = value.tolist()
        if isinstance(out, list):
            return out
        return [out]
    if scalar_missing(value):
        return []
    return [value]


def clean_str(value: Any, default: str = "") -> str:
    """Convert a scalar to a clean string, preserving empty default for missing."""
    if scalar_missing(value):
        return default
    return str(value)


def struct_get(value: Any, key: str, default: Any = None) -> Any:
    """Read a key from a struct-like value produced by pyarrow -> pandas."""
    if value is None:
        return default
    if isinstance(value, dict):
        return value.get(key, default)
    # Pandas may occasionally materialize structs as Series-like objects.
    try:
        return value[key]
    except Exception:
        return default


def nested_get(value: Any, keys: Sequence[str], default: Any = None) -> Any:
    """Read nested struct keys from a dict/Series-like value."""
    current = value
    for key in keys:
        current = struct_get(current, key, default=None)
        if current is None:
            return default
    return current


def collapse_list(value: Any, sep: str = "|") -> str:
    """Collapse a scalar/list-like value to a delimiter-separated string."""
    values = [clean_str(x) for x in as_list(value) if not scalar_missing(x) and clean_str(x) != ""]
    return sep.join(values)


# -----------------------------
# Gene list handling
# -----------------------------

def load_gene_list(path: Path) -> Dict[str, float]:
    """Load GeneList.txt as {HGNC symbol: globalScore}."""
    if not path.exists():
        raise FileNotFoundError(f"Gene list not found: {path}")
    df = pd.read_csv(path, sep="\t")
    required = {"symbol", "globalScore"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Gene list is missing required columns: {sorted(missing)}")
    out: Dict[str, float] = {}
    for _, row in df.iterrows():
        symbol = clean_str(row["symbol"]).strip()
        if not symbol:
            continue
        try:
            score = float(row["globalScore"])
        except Exception:
            score = float("nan")
        out[symbol] = score
    if not out:
        raise ValueError(f"Gene list contained no usable symbols: {path}")
    return out


def matched_genes_from_gencode(gencode_value: Any, gene_scores: Dict[str, float]) -> Tuple[List[str], List[float]]:
    """Return GeneList genes present in gencode.genes."""
    genes = [clean_str(x).strip() for x in as_list(struct_get(gencode_value, "genes", []))]
    genes = [g for g in genes if g]
    matched = [g for g in genes if g in gene_scores]
    scores = [gene_scores[g] for g in matched]
    return matched, scores


def best_gene(matched: Sequence[str], scores: Sequence[float]) -> Tuple[str, float]:
    """Pick the highest-scoring matched gene for deterministic single-gene fields."""
    if not matched:
        return "", float("nan")
    # If scores are NA, push them to the bottom but keep deterministic ordering.
    idx = max(range(len(matched)), key=lambda i: (-1 if math.isnan(scores[i]) else scores[i], matched[i]))
    return matched[idx], scores[idx]


# -----------------------------
# FAVOR-to-legacy annotation mapping
# -----------------------------

def map_legacy_exonic_category(consequence: Any, ref: Any = None, alt: Any = None) -> str:
    """
    Map FAVOR consequence-like strings to legacy IMPACT category terms.

    This is intentionally conservative. If no known mapping is found, returns the
    raw consequence string so we can inspect and expand the map later.
    """
    c = clean_str(consequence).strip()
    lc = c.lower().replace("_", " ").replace("-", " ")

    if not lc:
        return ""
    if any(x in lc for x in ["stopgain", "stop gained", "stop gained", "nonsense"]):
        return "stopgain"
    if any(x in lc for x in ["stoploss", "stop lost", "stop lost"]):
        return "stoploss"
    if "frameshift" in lc:
        # If ref/alt lengths are available, infer insertion/deletion direction.
        r = clean_str(ref)
        a = clean_str(alt)
        if r and a and len(a) > len(r):
            return "frameshift insertion"
        if r and a and len(r) > len(a):
            return "frameshift deletion"
        return "frameshift deletion"
    if any(x in lc for x in ["missense", "nonsynonymous"]):
        return "nonsynonymous SNV"
    if "inframe" in lc or "nonframeshift" in lc:
        r = clean_str(ref)
        a = clean_str(alt)
        if r and a and len(a) > len(r):
            return "nonframeshift insertion"
        if r and a and len(r) > len(a):
            return "nonframeshift deletion"
        return "nonframeshift insertion"
    return c


def gencode_info_string(matched_genes: Sequence[str]) -> str:
    """
    Create a legacy-compatible gene info string.

    Existing IMPACT prioritization code expects gene symbols and handles optional
    (dist=...) decorations, so we emit GENE(dist=0) for matched genes.
    """
    return ",".join(f"{gene}(dist=0)" for gene in matched_genes)


def flatten_annotation_row(row: pd.Series, gene_scores: Dict[str, float]) -> Optional[Dict[str, Any]]:
    """Flatten one FAVOR annotation row; return None if it has no GeneList match."""
    gencode = row.get("gencode")
    matched, scores = matched_genes_from_gencode(gencode, gene_scores)
    if not matched:
        return None

    best, best_score = best_gene(matched, scores)
    ref = row.get("ref_vcf")
    alt = row.get("alt_vcf")

    refseq = row.get("refseq")
    ucsc = row.get("ucsc")
    clinvar = row.get("clinvar")
    dbsnp = row.get("dbsnp")
    apc = row.get("apc")
    spliceai = row.get("spliceai")
    alphamissense = row.get("alphamissense")
    dbnsfp = row.get("dbnsfp")
    gnomad_genome = row.get("gnomad_genome")
    gnomad_exome = row.get("gnomad_exome")
    bravo = row.get("bravo")
    tg = row.get("tg")
    main = row.get("main")
    input_struct = row.get("input")

    gencode_consequence = struct_get(gencode, "consequence", "")
    refseq_consequence = struct_get(refseq, "consequence", "")
    ucsc_consequence = struct_get(ucsc, "consequence", "")

    clnsig = collapse_list(struct_get(clinvar, "clnsig", []), sep="\\")
    clndn = collapse_list(struct_get(clinvar, "clndn", []), sep="|")

    rsid = clean_str(struct_get(input_struct, "rsid", ""))
    if not rsid:
        rsid = clean_str(struct_get(dbsnp, "rsid", ""))

    out = {
        "chromosome": normalize_chrom(row.get("chromosome")),
        "position": int(row.get("position")),
        "ref": clean_str(ref),
        "alt": clean_str(alt),
        "allele": f"{clean_str(ref)},{clean_str(alt)}",
        "vid": row.get("vid"),
        "variant_vcf": clean_str(row.get("variant_vcf")),
        "rsid": rsid,
        "qual": clean_str(struct_get(input_struct, "qual", "")),
        "filter": clean_str(struct_get(input_struct, "filter", "")),
        "matched_gene": best,
        "matched_gene_score": best_score,
        "matched_gene_all": "|".join(matched),
        "matched_gene_score_all": "|".join("" if math.isnan(s) else str(s) for s in scores),
        "matched_gene_source": "gencode.genes",
        "gencode_genes": collapse_list(struct_get(gencode, "genes", []), sep="|"),
        "gencode_region_type": clean_str(struct_get(gencode, "region_type", "")),
        "gencode_consequence": clean_str(gencode_consequence),
        "refseq_consequence": clean_str(refseq_consequence),
        "ucsc_consequence": clean_str(ucsc_consequence),
        "genecode_comprehensive_info": gencode_info_string(matched),
        "genecode_comprehensive_exonic_category": map_legacy_exonic_category(gencode_consequence, ref, alt),
        "refseq_exonic_category": map_legacy_exonic_category(refseq_consequence, ref, alt),
        "ucsc_exonic_category": map_legacy_exonic_category(ucsc_consequence, ref, alt),
        "clnsig": clnsig,
        "clndn": clndn,
        "clinvar_gene": clean_str(struct_get(clinvar, "gene", "")),
        "apc_protein_function_v3": nested_get(apc, ["protein_function_v3"], float("nan")),
        "spliceai_max_ds": nested_get(spliceai, ["max_ds"], float("nan")),
        "alphamissense_max_pathogenicity": nested_get(alphamissense, ["max_pathogenicity"], float("nan")),
        "revel": nested_get(dbnsfp, ["revel"], float("nan")),
        "gnomad_genome_af": nested_get(gnomad_genome, ["af"], float("nan")),
        "gnomad_exome_af": nested_get(gnomad_exome, ["af"], float("nan")),
        "bravo_af": nested_get(bravo, ["bravo_af"], float("nan")),
        "tg_all": nested_get(tg, ["tg_all"], float("nan")),
        "cadd_phred": nested_get(main, ["cadd", "phred"], float("nan")),
    }
    return out


# -----------------------------
# Parquet loading
# -----------------------------

def parquet_path(base_dir: Path, chrom: str) -> Path:
    """Return chromosome parquet path."""
    chrom = normalize_chrom(chrom)
    path = base_dir / f"chromosome={chrom}" / "data.parquet"
    if not path.exists():
        raise FileNotFoundError(f"Missing chromosome parquet: {path}")
    return path


def read_parquet_file(path: Path, columns: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """
    Read one physical parquet file without Hive partition inference.

    This avoids pyarrow dataset partition type conflicts when the path includes
    chromosome=... and the file also contains a chromosome column.
    """
    return pq.ParquetFile(path).read(columns=list(columns) if columns else None).to_pandas()


def load_samples(genotypes_dir: Path) -> List[str]:
    """Load sample order from genotypes_dir/samples.txt."""
    samples_path = genotypes_dir / "samples.txt"
    if not samples_path.exists():
        raise FileNotFoundError(f"Missing samples.txt: {samples_path}")
    samples = [line.strip() for line in samples_path.read_text().splitlines() if line.strip()]
    if not samples:
        raise ValueError(f"No sample IDs found in {samples_path}")
    return samples


def dosage_at(value: Any, sample_index: int) -> float:
    """Extract one sample dosage from a fixed-size list value."""
    values = as_list(value)
    if sample_index >= len(values):
        raise IndexError(f"sample_index={sample_index} but dosage vector length is {len(values)}")
    x = values[sample_index]
    if scalar_missing(x):
        return float("nan")
    return float(x)


# -----------------------------
# Main flattening logic
# -----------------------------

def flatten_one(
    annotated_dir: Path,
    genotypes_dir: Path,
    gene_list: Path,
    sample_id: str,
    chromosome: str,
    dosage_threshold: float,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Flatten one sample/chromosome into a filtered table plus summary."""
    chrom = normalize_chrom(chromosome)
    gene_scores = load_gene_list(gene_list)
    samples = load_samples(genotypes_dir)
    if sample_id not in samples:
        raise ValueError(f"sample_id {sample_id!r} not found in samples.txt: {samples}")
    sample_index = samples.index(sample_id)

    geno_path = parquet_path(genotypes_dir, chrom)
    ann_path = parquet_path(annotated_dir, chrom)

    geno_cols = ["chromosome", "position", "ref", "alt", "maf", "dosages"]
    geno = read_parquet_file(geno_path, columns=geno_cols)
    total_genotype_rows = len(geno)
    geno["chromosome"] = geno["chromosome"].map(normalize_chrom)
    geno["dosage"] = geno["dosages"].map(lambda x: dosage_at(x, sample_index))
    carried = geno[geno["dosage"] > dosage_threshold].copy()
    carried_count = len(carried)
    carried = carried.drop(columns=["dosages"])

    ann_cols = [
        "input", "vid", "chromosome", "position", "ref_vcf", "alt_vcf", "variant_vcf",
        "gencode", "ucsc", "refseq", "dbsnp", "clinvar", "apc", "spliceai",
        "alphamissense", "dbnsfp", "gnomad_genome", "gnomad_exome", "bravo", "tg", "main",
    ]
    ann = read_parquet_file(ann_path, columns=ann_cols)
    total_annotation_rows = len(ann)

    flattened_rows: List[Dict[str, Any]] = []
    for _, row in ann.iterrows():
        flat = flatten_annotation_row(row, gene_scores)
        if flat is not None:
            flattened_rows.append(flat)
    ann_flat = pd.DataFrame(flattened_rows)
    gene_matching_annotation_rows = len(ann_flat)

    if ann_flat.empty or carried.empty:
        out = pd.DataFrame()
    else:
        carried_key = carried.rename(columns={"ref": "ref", "alt": "alt"})
        ann_flat["chromosome"] = ann_flat["chromosome"].map(normalize_chrom)
        out = carried_key.merge(
            ann_flat,
            on=["chromosome", "position", "ref", "alt"],
            how="inner",
            suffixes=("_geno", ""),
        )

    if not out.empty:
        out.insert(0, "sample_id", sample_id)
        out.insert(1, "variant_id", range(1, len(out) + 1))
        # Keep preferred column order for easier validation and GDS writing.
        preferred = [
            "sample_id", "variant_id", "chromosome", "position", "ref", "alt", "allele",
            "dosage", "maf", "vid", "variant_vcf", "rsid", "qual", "filter",
            "matched_gene", "matched_gene_score", "matched_gene_all", "matched_gene_score_all", "matched_gene_source",
            "gencode_genes", "gencode_region_type", "gencode_consequence",
            "genecode_comprehensive_info", "genecode_comprehensive_exonic_category",
            "refseq_exonic_category", "ucsc_exonic_category",
            "clnsig", "clndn", "clinvar_gene", "apc_protein_function_v3",
            "spliceai_max_ds", "alphamissense_max_pathogenicity", "revel",
            "gnomad_genome_af", "gnomad_exome_af", "bravo_af", "tg_all", "cadd_phred",
        ]
        remaining = [c for c in out.columns if c not in preferred]
        out = out[[c for c in preferred if c in out.columns] + remaining]

    summary = {
        "sample_id": sample_id,
        "sample_index": sample_index,
        "samples": samples,
        "chromosome": chrom,
        "gene_list": str(gene_list),
        "gene_count": len(gene_scores),
        "dosage_threshold": dosage_threshold,
        "genotype_path": str(geno_path),
        "annotation_path": str(ann_path),
        "total_genotype_rows": total_genotype_rows,
        "carried_variants_for_sample": carried_count,
        "total_annotation_rows": total_annotation_rows,
        "gene_matching_annotation_rows": gene_matching_annotation_rows,
        "final_carried_gene_matched_rows": len(out),
        "rows_with_clnsig": int((out["clnsig"].fillna("") != "").sum()) if not out.empty and "clnsig" in out else 0,
        "rows_with_apc_protein_function_v3": int(out["apc_protein_function_v3"].notna().sum()) if not out.empty and "apc_protein_function_v3" in out else 0,
        "matched_genes_observed": sorted(out["matched_gene"].dropna().unique().tolist()) if not out.empty and "matched_gene" in out else [],
    }
    return out, summary


# -----------------------------
# CLI
# -----------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a flat preview/intermediate table from FAVOR annotated + genotype parquet for one sample/chromosome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--annotated-dir", required=True, type=Path, help="FAVOR annotated output directory")
    parser.add_argument("--genotypes-dir", required=True, type=Path, help="FAVOR genotype output directory")
    parser.add_argument("--gene-list", required=True, type=Path, help="GeneList.txt with symbol and globalScore columns")
    parser.add_argument("--sample-id", required=True, help="Sample ID to extract, as listed in genotypes_dir/samples.txt")
    parser.add_argument("--chromosome", required=True, help="Chromosome to process, e.g. 1, 22, X, Y")
    parser.add_argument("--dosage-threshold", type=float, default=0.0, help="Keep variants with dosage > threshold")
    parser.add_argument("--out", required=True, type=Path, help="Output flattened parquet path")
    parser.add_argument("--preview-csv", type=Path, help="Optional CSV containing the first --preview-rows rows")
    parser.add_argument("--preview-rows", type=int, default=50, help="Rows to write to --preview-csv")
    parser.add_argument("--summary-json", type=Path, help="Optional JSON summary path")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        flat, summary = flatten_one(
            annotated_dir=args.annotated_dir,
            genotypes_dir=args.genotypes_dir,
            gene_list=args.gene_list,
            sample_id=args.sample_id,
            chromosome=args.chromosome,
            dosage_threshold=args.dosage_threshold,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    flat.to_parquet(args.out, index=False, engine="pyarrow")

    if args.preview_csv:
        args.preview_csv.parent.mkdir(parents=True, exist_ok=True)
        flat.head(args.preview_rows).to_csv(args.preview_csv, index=False)

    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    print("Flatten preview completed")
    print(f"  output:                          {args.out}")
    print(f"  sample_id:                       {summary['sample_id']}")
    print(f"  chromosome:                      {summary['chromosome']}")
    print(f"  gene_count:                      {summary['gene_count']}")
    print(f"  total_genotype_rows:             {summary['total_genotype_rows']}")
    print(f"  carried_variants_for_sample:     {summary['carried_variants_for_sample']}")
    print(f"  total_annotation_rows:           {summary['total_annotation_rows']}")
    print(f"  gene_matching_annotation_rows:   {summary['gene_matching_annotation_rows']}")
    print(f"  final_carried_gene_matched_rows: {summary['final_carried_gene_matched_rows']}")
    print(f"  rows_with_clnsig:                {summary['rows_with_clnsig']}")
    print(f"  rows_with_apc_protein_function:  {summary['rows_with_apc_protein_function_v3']}")
    if summary["matched_genes_observed"]:
        preview = ", ".join(summary["matched_genes_observed"][:20])
        if len(summary["matched_genes_observed"]) > 20:
            preview += ", ..."
        print(f"  matched_genes_observed:          {preview}")
    else:
        print("  matched_genes_observed:          none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
