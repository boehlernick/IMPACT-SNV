#!/usr/bin/env python3
"""Flatten FAVOR annotation + IMPACT-SNV genotype parquet for GDS building.

This script creates one compact, per-sample, per-chromosome flat parquet table
from:

    annotated_dir/chromosome=<chrom>/data.parquet
    genotypes_dir/chromosome=<chrom>/data.parquet
    genotypes_dir/samples.txt
    GeneList.txt

The output is consumed by impact_snv/resources/favor_flat_to_seqarray_gds.R.
It intentionally emits both modern IMPACT-SNV fields and the legacy field names
required by the current R GDS writer.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd
import pyarrow.parquet as pq

from impact_snv.gds import contract


def normalize_chrom(value: Any) -> str:
    if value is None:
        return ""
    text = str(value)
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.endswith(".0"):
        text = text[:-2]
    return text


def scalar_missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except Exception:
        return False


def as_list(value: Any) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if hasattr(value, "tolist"):
        out = value.tolist()
        return out if isinstance(out, list) else [out]
    if scalar_missing(value):
        return []
    return [value]


def clean_str(value: Any, default: str = "") -> str:
    if scalar_missing(value):
        return default
    return str(value)


def numeric_or_none(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    try:
        out = float(value)
    except Exception:
        return None
    if math.isnan(out):
        return None
    return out


def numeric_or_nan(value: Any) -> float:
    out = numeric_or_none(value)
    return float("nan") if out is None else float(out)


def numeric_or_zero(value: Any) -> float:
    out = numeric_or_none(value)
    return 0.0 if out is None else float(out)


def allele_frequency_to_maf(value: Any) -> Optional[float]:
    af = numeric_or_none(value)
    if af is None or af < 0 or af > 1:
        return None
    return float(min(af, 1.0 - af))


def struct_get(value: Any, key: str, default: Any = None) -> Any:
    if value is None:
        return default
    if isinstance(value, dict):
        return value.get(key, default)
    try:
        return value[key]
    except Exception:
        return default


def nested_get(value: Any, keys: Sequence[str], default: Any = None) -> Any:
    current = value
    for key in keys:
        current = struct_get(current, key, default=None)
        if current is None:
            return default
    return current


def collapse_list(value: Any, sep: str = ";") -> str:
    vals = []
    for x in as_list(value):
        sx = clean_str(x)
        if sx:
            vals.append(sx)
    return sep.join(vals)


def normalize_exonic_category_chunks(value: Any) -> List[str]:
    text = collapse_list(value, sep=";").lower()
    if not text:
        return []
    chunks: List[str] = []
    for chunk in re.split(r"[|,;/]+", text):
        normalized = re.sub(r"[_\-]+", " ", chunk)
        normalized = re.sub(r"\s+", " ", normalized).strip()
        if normalized:
            chunks.append(normalized)
    return chunks


def infer_indel_direction(ref: Any, alt: Any, text: str = "") -> str:
    if re.search(r"\binsertion\b", text):
        return "insertion"
    if re.search(r"\bdeletion\b", text):
        return "deletion"
    ref_text = clean_str(ref).strip()
    alt_text = clean_str(alt).strip()
    if not ref_text or not alt_text or "," in alt_text:
        return ""
    if len(alt_text) > len(ref_text):
        return "insertion"
    if len(alt_text) < len(ref_text):
        return "deletion"
    return ""


def classify_exonic_category_labels(consequence: Any, ref: Any = None, alt: Any = None) -> List[str]:
    labels: List[str] = []
    for chunk in normalize_exonic_category_chunks(consequence):
        direction = infer_indel_direction(ref, alt, chunk)
        if re.search(r"\bstop\s*gain(?:ed)?\b|\bnonsense\b", chunk):
            labels.append("stopgain")
        if re.search(r"\bstop\s*loss\b|\bstop\s*lost\b", chunk):
            labels.append("stoploss")
        if re.search(r"\bframeshift\b", chunk):
            labels.append(f"frameshift {direction}".strip())
        if re.search(r"\bnonframeshift\b|\binframe\b", chunk):
            labels.append(f"nonframeshift {direction}".strip())
        if re.search(r"\bnonsynonymous\b|\bmissense\b", chunk):
            labels.append("nonsynonymous SNV")
        if re.search(r"\bsynonymous\b", chunk):
            labels.append("synonymous SNV")
        if re.search(r"\bsplic", chunk):
            labels.append("splicing")
    seen = set()
    ordered: List[str] = []
    for label in labels:
        if label and label not in seen:
            seen.add(label)
            ordered.append(label)
    return ordered


def load_gene_list(path: Path) -> Dict[str, float]:
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
        out[symbol] = numeric_or_nan(row["globalScore"])
    if not out:
        raise ValueError(f"Gene list contained no usable symbols: {path}")
    return out


def matched_genes_from_gencode(gencode_value: Any, gene_scores: Dict[str, float]) -> Tuple[List[str], List[float]]:
    genes = [clean_str(x).strip() for x in as_list(struct_get(gencode_value, "genes", []))]
    genes = [g for g in genes if g]
    matched = [g for g in genes if g in gene_scores]
    scores = [gene_scores[g] for g in matched]
    return matched, scores


def best_gene(matched: Sequence[str], scores: Sequence[float]) -> Tuple[str, float]:
    if not matched:
        return "", float("nan")
    idx = max(range(len(matched)), key=lambda i: (-1 if math.isnan(scores[i]) else scores[i], matched[i]))
    return matched[idx], float(scores[idx])


def extract_maf_from_annotation_row(row: pd.Series) -> Tuple[float, bool]:
    candidates = [
        row.get("bravo_af"),
        row.get("gnomad_genome_af"),
        row.get("gnomad_exome_af"),
        nested_get(row.get("bravo"), ["bravo_af"]),
        nested_get(row.get("gnomad_genome"), ["af"]),
        nested_get(row.get("gnomad_exome"), ["af"]),
    ]
    for value in candidates:
        maf = allele_frequency_to_maf(value)
        if maf is not None:
            return maf, False
    return 0.0, True


def extract_af(row: pd.Series, top_level: str, nested_col: str, nested_key: str) -> float:
    value = row.get(top_level)
    if value is None:
        value = nested_get(row.get(nested_col), [nested_key])
    return numeric_or_nan(value)


def first_nonempty(*values: Any, default: str = "") -> str:
    for value in values:
        s = clean_str(value)
        if s:
            return s
    return default


def map_legacy_exonic_category(consequence: Any, ref: Any = None, alt: Any = None) -> str:
    text = collapse_list(consequence, sep=";").lower()
    if not text:
        return "unknown"
    labels = classify_exonic_category_labels(consequence, ref, alt)
    for label in (
        "stopgain",
        "frameshift deletion",
        "frameshift insertion",
        "frameshift",
        "stoploss",
        "nonsynonymous SNV",
        "synonymous SNV",
        "splicing",
        "nonframeshift deletion",
        "nonframeshift insertion",
        "nonframeshift",
    ):
        if label in labels:
            return label
    return text[:200]


def get_variant_fields(row: pd.Series) -> Tuple[str, int, str, str]:
    chrom = normalize_chrom(row.get("chromosome", row.get("chrom", "")))
    pos = numeric_or_none(row.get("position"))
    if pos is None:
        pos = numeric_or_none(row.get("pos"))
    ref = clean_str(row.get("ref_vcf", row.get("ref", "")))
    alt = clean_str(row.get("alt_vcf", row.get("alt", "")))
    return chrom, int(pos) if pos is not None else -1, ref, alt


def extract_clnsig(row: pd.Series) -> str:
    clinvar = row.get("clinvar")
    return first_nonempty(
        row.get("clnsig"),
        row.get("clinvar_clnsig"),
        nested_get(clinvar, ["clinical_significance"]),
        nested_get(clinvar, ["clnsig"]),
        nested_get(clinvar, ["CLNSIG"]),
    )


def extract_apc_protein_function(row: pd.Series) -> Tuple[float, bool]:
    apc = row.get("apc")
    candidates = [
        row.get("apc_protein_function_v3"),
        row.get("apc_protein_function"),
        nested_get(apc, ["protein_function_v3"]),
        nested_get(apc, ["protein_function"]),
        nested_get(apc, ["score"]),
        nested_get(apc, ["v3"]),
    ]
    for value in candidates:
        out = numeric_or_none(value)
        if out is not None:
            return float(out), False
    return 0.0, True


def flatten_annotation_row(row: pd.Series, gene_scores: Dict[str, float]) -> Optional[Dict[str, Any]]:
    gencode = row.get("gencode")
    matched, scores = matched_genes_from_gencode(gencode, gene_scores)
    if not matched:
        return None

    chrom, pos, ref, alt = get_variant_fields(row)
    if pos < 0 or not ref or not alt:
        return None

    best, best_score = best_gene(matched, scores)
    consequence = struct_get(gencode, "consequence")
    region_type = clean_str(struct_get(gencode, "region_type"))
    exonic_category = map_legacy_exonic_category(consequence, ref, alt)
    genes_joined = ";".join(matched)
    scores_joined = ";".join("" if math.isnan(s) else str(float(s)) for s in scores)

    bravo_af = extract_af(row, "bravo_af", "bravo", "bravo_af")
    gnomad_genome_af = extract_af(row, "gnomad_genome_af", "gnomad_genome", "af")
    gnomad_exome_af = extract_af(row, "gnomad_exome_af", "gnomad_exome", "af")
    vid = first_nonempty(row.get("vid"), default=f"{chrom}:{pos}:{ref}:{alt}")
    variant_vcf = first_nonempty(row.get("variant_vcf"), default=f"{chrom}:{pos}:{ref}:{alt}")
    transcript_info = collapse_list(struct_get(gencode, "transcripts"), sep=";")
    clnsig = extract_clnsig(row)
    apc_pf, apc_fallback = extract_apc_protein_function(row)
    maf, maf_fallback = extract_maf_from_annotation_row(row)

    fallback_flags: List[str] = []
    if not transcript_info:
        fallback_flags.append("GENECODE_INFO_FALLBACK_TO_MATCHED_GENES")
    if not clnsig:
        fallback_flags.append("CLNSIG_MISSING_DEFAULT_EMPTY")
    if apc_fallback:
        fallback_flags.append("APC_PROTEIN_FUNCTION_DEFAULT_0")
    if maf_fallback:
        fallback_flags.append("MAF_DEFAULT_0")

    return {
        "chromosome": chrom,
        "position": pos,
        "ref": ref,
        "alt": alt,
        "allele": alt,
        # variant_id is filled with a sequential integer after sample/chrom filtering.
        "vid": vid,
        "variant_vcf": variant_vcf,
        "maf": maf,
        "bravo_af": bravo_af,
        "gnomad_genome_af": gnomad_genome_af,
        "gnomad_exome_af": gnomad_exome_af,
        "gene": best,
        "genes": genes_joined,
        "matched_genes": genes_joined,
        "matched_gene": best,
        "matched_gene_score": best_score,
        "matched_gene_all": genes_joined,
        "matched_gene_score_all": scores_joined,
        "matched_gene_source": "gencode.genes",
        "globalScore": best_score,
        "global_score": best_score,
        "gencode_genes": genes_joined,
        "gencode_region_type": region_type,
        "gencode_consequence": collapse_list(consequence, sep=";"),
        "genecode_comprehensive_info": transcript_info if transcript_info else genes_joined,
        "genecode_comprehensive_exonic_category": exonic_category,
        "refseq_exonic_category": exonic_category,
        "ucsc_exonic_category": exonic_category,
        "Func_refGene": region_type,
        "Gene_refGene": best,
        "GeneDetail_refGene": genes_joined,
        "ExonicFunc_refGene": exonic_category,
        "AAChange_refGene": transcript_info,
        "clnsig": clnsig,
        "clinvar_clnsig": clnsig,
        "apc_protein_function": apc_pf,
        "apc_protein_function_v3": apc_pf,
        "cadd_phred": numeric_or_nan(nested_get(row.get("main"), ["cadd", "phred"])),
        "impact_fallback_flags": ";".join(fallback_flags),
        "impact_fallback_count": len(fallback_flags),
    }


def parquet_path(base_dir: Path, chrom: str) -> Path:
    chrom = normalize_chrom(chrom)
    path = base_dir / f"chromosome={chrom}" / "data.parquet"
    if not path.exists():
        raise FileNotFoundError(f"Missing chromosome parquet: {path}")
    return path


def read_parquet_file(path: Path, columns: Optional[Sequence[str]] = None) -> pd.DataFrame:
    # Important: ParquetFile.read() avoids Hive partition inference from chromosome=*/.
    pf = pq.ParquetFile(path)
    table = pf.read(columns=list(columns) if columns else None)
    return table.to_pandas()


def load_samples(genotypes_dir: Path) -> List[str]:
    samples_path = genotypes_dir / "samples.txt"
    if not samples_path.exists():
        raise FileNotFoundError(f"Missing samples.txt: {samples_path}")
    samples = [line.strip() for line in samples_path.read_text().splitlines() if line.strip()]
    if not samples:
        raise ValueError(f"No sample IDs found in {samples_path}")
    return samples


def dosage_at(value: Any, sample_index: int) -> float:
    values = as_list(value)
    if sample_index >= len(values):
        raise IndexError(f"sample_index={sample_index} but dosage vector length is {len(values)}")
    x = values[sample_index]
    if scalar_missing(x):
        return float("nan")
    return float(x)


def normalize_key_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "chromosome" not in out.columns:
        if "chrom" in out.columns:
            out["chromosome"] = out["chrom"]
        elif "CHROM" in out.columns:
            out["chromosome"] = out["CHROM"]
    if "position" not in out.columns:
        if "pos" in out.columns:
            out["position"] = out["pos"]
        elif "POS" in out.columns:
            out["position"] = out["POS"]
    if "ref" not in out.columns:
        if "ref_vcf" in out.columns:
            out["ref"] = out["ref_vcf"]
        elif "REF" in out.columns:
            out["ref"] = out["REF"]
    if "alt" not in out.columns:
        if "alt_vcf" in out.columns:
            out["alt"] = out["alt_vcf"]
        elif "ALT" in out.columns:
            out["alt"] = out["ALT"]
    for col in ("chromosome", "position", "ref", "alt"):
        if col not in out.columns:
            raise KeyError(col)
    out["chromosome"] = out["chromosome"].map(normalize_chrom)
    out["position"] = out["position"].astype("int64")
    out["ref"] = out["ref"].map(clean_str)
    out["alt"] = out["alt"].map(clean_str)
    return out


def flatten_one(
    annotated_dir: Path,
    genotypes_dir: Path,
    gene_list: Path,
    sample_id: str,
    chromosome: str,
    dosage_threshold: float,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    chrom = normalize_chrom(chromosome)
    gene_scores = load_gene_list(gene_list)
    samples = load_samples(genotypes_dir)
    if sample_id not in samples:
        raise ValueError(f"sample_id {sample_id!r} not found in samples.txt: {samples}")
    sample_index = samples.index(sample_id)

    ann_raw = read_parquet_file(parquet_path(annotated_dir, chrom))
    geno_raw = read_parquet_file(parquet_path(genotypes_dir, chrom))

    flat_ann_rows = []
    for _, row in ann_raw.iterrows():
        rec = flatten_annotation_row(row, gene_scores)
        if rec is not None:
            flat_ann_rows.append(rec)
    ann = pd.DataFrame(flat_ann_rows) if flat_ann_rows else pd.DataFrame()
    if not ann.empty:
        ann = normalize_key_columns(ann)

    geno = normalize_key_columns(geno_raw)
    dosage_col = "dosages" if "dosages" in geno.columns else "dosage"
    if dosage_col not in geno.columns:
        raise KeyError("dosages")
    geno["dosage"] = geno[dosage_col].map(lambda x: dosage_at(x, sample_index))
    geno = geno[~geno["dosage"].isna()].copy()
    carried = geno[geno["dosage"] > dosage_threshold].copy()

    join_cols = ["chromosome", "position", "ref", "alt"]
    if ann.empty or carried.empty:
        merged = pd.DataFrame(columns=required_output_columns())
    else:
        merged = carried.merge(ann, on=join_cols, how="inner", suffixes=("_gt", ""))

    merged = finalize_flat_schema(merged, sample_id)

    matched_genes = sorted(set(g for value in ann.get("matched_genes", []) for g in clean_str(value).split(";") if g)) if not ann.empty else []
    rows_with_clnsig = int((merged.get("clnsig", pd.Series(dtype=str)).fillna("").astype(str) != "").sum()) if len(merged) else 0
    rows_with_apc = int((merged.get("apc_protein_function_v3", pd.Series(dtype=float)).fillna(0).astype(float) != 0).sum()) if len(merged) else 0
    fallback_series = merged.get("impact_fallback_flags", pd.Series(dtype=str)).fillna("").astype(str) if len(merged) else pd.Series(dtype=str)
    rows_with_fallbacks = int((fallback_series != "").sum()) if len(merged) else 0
    fallback_counts = Counter(
        flag
        for value in fallback_series
        for flag in value.split(";")
        if flag
    )
    summary = {
        "sample_id": sample_id,
        "chromosome": chrom,
        "gene_count": len(gene_scores),
        "total_genotype_rows": int(len(geno_raw)),
        "carried_variants_for_sample": int(len(carried)),
        "total_annotation_rows": int(len(ann_raw)),
        "gene_matching_annotation_rows": int(len(ann)),
        "final_carried_gene_matched_rows": int(len(merged)),
        "rows_with_clnsig": rows_with_clnsig,
        "rows_with_apc_protein_function": rows_with_apc,
        "rows_with_compatibility_fallbacks": rows_with_fallbacks,
        "compatibility_fallback_counts": dict(sorted(fallback_counts.items())),
        "matched_genes_observed": matched_genes[:200],
    }
    return merged, summary


def required_output_columns() -> List[str]:
    """Return the canonical required flat columns from the package contract."""
    return contract.required_output_columns()


def finalize_flat_schema(df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    out = df.copy()
    n = len(out)

    if "impact_fallback_flags" not in out.columns:
        out["impact_fallback_flags"] = ""

    def append_fallback_flag(flag: str) -> None:
        if len(out) == 0:
            return
        values = out["impact_fallback_flags"].fillna("").astype(str)
        out["impact_fallback_flags"] = values.map(lambda x: flag if not x else f"{x};{flag}")

    if "sample_id" not in out.columns:
        out["sample_id"] = sample_id
    if "allele" not in out.columns:
        out["allele"] = out["alt"] if "alt" in out.columns else ""
        append_fallback_flag("BACKFILL_ALLELE")
    if "maf" not in out.columns:
        out["maf"] = 0.0
        append_fallback_flag("BACKFILL_MAF_DEFAULT_0")
    if "vid" not in out.columns:
        out["vid"] = ""
        append_fallback_flag("BACKFILL_VID_EMPTY")
    if "variant_vcf" not in out.columns:
        if all(c in out.columns for c in ["chromosome", "position", "ref", "alt"]):
            out["variant_vcf"] = out["chromosome"].astype(str) + ":" + out["position"].astype(str) + ":" + out["ref"].astype(str) + ":" + out["alt"].astype(str)
            append_fallback_flag("BACKFILL_VARIANT_VCF_FROM_KEYS")
        else:
            out["variant_vcf"] = ""
            append_fallback_flag("BACKFILL_VARIANT_VCF_EMPTY")
    if "matched_gene" not in out.columns:
        out["matched_gene"] = out["gene"] if "gene" in out.columns else ""
        append_fallback_flag("BACKFILL_MATCHED_GENE")
    if "matched_gene_score" not in out.columns:
        out["matched_gene_score"] = out["globalScore"] if "globalScore" in out.columns else float("nan")
        append_fallback_flag("BACKFILL_MATCHED_GENE_SCORE")
    if "genecode_comprehensive_info" not in out.columns:
        out["genecode_comprehensive_info"] = out["GeneDetail_refGene"] if "GeneDetail_refGene" in out.columns else ""
        append_fallback_flag("BACKFILL_GENECODE_INFO")
    if "genecode_comprehensive_exonic_category" not in out.columns:
        out["genecode_comprehensive_exonic_category"] = out["ExonicFunc_refGene"] if "ExonicFunc_refGene" in out.columns else "unknown"
        append_fallback_flag("BACKFILL_GENECODE_EXONIC_CATEGORY")
    if "refseq_exonic_category" not in out.columns:
        out["refseq_exonic_category"] = out["genecode_comprehensive_exonic_category"]
        append_fallback_flag("BACKFILL_REFSEQ_EXONIC_CATEGORY")
    if "ucsc_exonic_category" not in out.columns:
        out["ucsc_exonic_category"] = out["genecode_comprehensive_exonic_category"]
        append_fallback_flag("BACKFILL_UCSC_EXONIC_CATEGORY")
    if "clnsig" not in out.columns:
        out["clnsig"] = ""
        append_fallback_flag("BACKFILL_CLNSIG_EMPTY")
    if "apc_protein_function_v3" not in out.columns:
        out["apc_protein_function_v3"] = 0.0
        append_fallback_flag("BACKFILL_APC_DEFAULT_0")

    # variant_id must be integer-like for the current R writer. Assign sequential
    # IDs per flat file; the R writer may reassign global IDs after read.
    out["variant_id"] = range(1, n + 1)

    # Normalize required scalar types.
    str_cols = ["sample_id", "chromosome", "ref", "alt", "allele", "vid", "variant_vcf", "matched_gene",
                "genecode_comprehensive_info", "genecode_comprehensive_exonic_category", "refseq_exonic_category",
                "ucsc_exonic_category", "clnsig"]
    for col in str_cols:
        if col in out.columns:
            out[col] = out[col].fillna("").astype(str)
    out["impact_fallback_flags"] = out["impact_fallback_flags"].fillna("").astype(str)
    out["impact_fallback_flags"] = out["impact_fallback_flags"].map(
        lambda x: ";".join(dict.fromkeys([flag for flag in x.split(";") if flag]))
    )
    out["impact_fallback_count"] = out["impact_fallback_flags"].map(lambda x: 0 if not x else len(x.split(";"))).astype("int64")
    for col in ["position", "variant_id"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0).astype("int64")
    for col in ["dosage", "maf", "matched_gene_score", "apc_protein_function_v3"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
            if col in ["dosage", "maf", "apc_protein_function_v3"]:
                out[col] = out[col].fillna(0.0)

    first = required_output_columns()
    optional_first = [
        "impact_fallback_flags", "impact_fallback_count",
        "matched_gene_all", "matched_gene_score_all", "matched_gene_source", "gene", "genes", "matched_genes",
        "globalScore", "global_score", "Func_refGene", "Gene_refGene", "GeneDetail_refGene", "ExonicFunc_refGene",
        "AAChange_refGene", "clinvar_clnsig", "apc_protein_function", "bravo_af", "gnomad_genome_af",
        "gnomad_exome_af", "tg_all", "gencode_genes", "gencode_region_type", "gencode_consequence", "cadd_phred",
    ]
    ordered = [c for c in first + optional_first if c in out.columns] + [c for c in out.columns if c not in first + optional_first]
    return out[ordered]


def write_outputs(df: pd.DataFrame, summary: Dict[str, Any], out: Path, preview_csv: Optional[Path], preview_rows: int, summary_json: Optional[Path]) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    # Validate schema against the package contract before writing.
    contract.validate_flat_df(df)
    df.to_parquet(out, index=False, engine="pyarrow")
    if preview_csv:
        preview_csv.parent.mkdir(parents=True, exist_ok=True)
        df.head(preview_rows).to_csv(preview_csv, index=False)
    if summary_json:
        summary_json.parent.mkdir(parents=True, exist_ok=True)
        summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")


def print_summary(out: Path, summary: Dict[str, Any]) -> None:
    print("Flatten preview completed")
    print(f"output:                          {out}")
    for key in [
        "sample_id", "chromosome", "gene_count", "total_genotype_rows", "carried_variants_for_sample",
        "total_annotation_rows", "gene_matching_annotation_rows", "final_carried_gene_matched_rows",
        "rows_with_clnsig", "rows_with_apc_protein_function",
    ]:
        print(f"{key + ':':32} {summary.get(key)}")
    genes = summary.get("matched_genes_observed") or []
    text = ", ".join(genes[:20]) + (", ..." if len(genes) > 20 else "") if genes else "none"
    print(f"matched_genes_observed:          {text}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a flat preview/intermediate table from FAVOR annotated + genotype parquet for one sample/chromosome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--annotated-dir", required=True, type=Path)
    parser.add_argument("--genotypes-dir", required=True, type=Path)
    parser.add_argument("--gene-list", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--dosage-threshold", type=float, default=0.0)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--preview-csv", type=Path)
    parser.add_argument("--preview-rows", type=int, default=50)
    parser.add_argument("--summary-json", type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        df, summary = flatten_one(
            annotated_dir=args.annotated_dir,
            genotypes_dir=args.genotypes_dir,
            gene_list=args.gene_list,
            sample_id=args.sample_id,
            chromosome=args.chromosome,
            dosage_threshold=args.dosage_threshold,
        )
        write_outputs(df, summary, args.out, args.preview_csv, args.preview_rows, args.summary_json)
        print_summary(args.out, summary)
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
