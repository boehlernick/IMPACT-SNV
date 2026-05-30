"""GDS flat file contract for IMPACT-SNV.

This module centralizes the required columns and optional fields expected by
the flat parquet produced by `impact_snv.gds.flatten` and consumed by the R
GDS writer `impact_snv.resources.favor_flat_to_seqarray_gds.R`.

Do not change these lists silently — update tests and R writer if contract
changes are required.
"""

from __future__ import annotations

from typing import List
from typing import Iterable


REQUIRED_FLAT_COLUMNS: List[str] = [
    "sample_id",
    "variant_id",
    "chromosome",
    "position",
    "ref",
    "alt",
    "allele",
    "dosage",
    "maf",
    "vid",
    "variant_vcf",
    "matched_gene",
    "matched_gene_score",
    "genecode_comprehensive_info",
    "genecode_comprehensive_exonic_category",
    "refseq_exonic_category",
    "ucsc_exonic_category",
    "clnsig",
    "apc_protein_function_v3",
]


OPTIONAL_FLAT_COLUMNS: List[str] = [
    "impact_fallback_flags",
    "impact_fallback_count",
    "matched_gene_all",
    "matched_gene_score_all",
    "matched_gene_source",
    "gene",
    "genes",
    "matched_genes",
    "globalScore",
    "global_score",
    "Func_refGene",
    "Gene_refGene",
    "GeneDetail_refGene",
    "ExonicFunc_refGene",
    "AAChange_refGene",
    "clinvar_clnsig",
    "apc_protein_function",
    "bravo_af",
    "gnomad_genome_af",
    "gnomad_exome_af",
    "tg_all",
    "gencode_genes",
    "gencode_region_type",
    "gencode_consequence",
    "cadd_phred",
]


def required_output_columns() -> List[str]:
    """Return the canonical ordered list of required flat columns."""
    return list(REQUIRED_FLAT_COLUMNS)


def validate_flat_df_columns(columns: Iterable[str]) -> List[str]:
    """Return list of missing required columns for the provided column iterable.

    Parameters
    - columns: Iterable of column names (e.g., `df.columns`).

    Returns
    - list of missing required column names (empty if none missing).
    """
    colset = set(columns)
    missing = [c for c in REQUIRED_FLAT_COLUMNS if c not in colset]
    return missing


def validate_flat_df(df, raise_on_missing: bool = True) -> List[str]:
    """Validate a flat dataframe against the contract.

    If `raise_on_missing` is True a ``ValueError`` will be raised when required
    columns are missing; otherwise the list of missing columns is returned.
    """
    missing = validate_flat_df_columns(getattr(df, "columns", []))
    if missing and raise_on_missing:
        raise ValueError(f"Flat dataframe missing required columns: {missing}")
    return missing
