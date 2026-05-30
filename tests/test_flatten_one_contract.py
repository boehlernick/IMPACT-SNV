from pathlib import Path

import pandas as pd
import pytest

from impact_snv.gds import contract
from impact_snv.gds.flatten import flatten_one, map_legacy_exonic_category


def _write_fixture_parquet(base: Path, chrom: str, df: pd.DataFrame) -> None:
    out_dir = base / f"chromosome={chrom}"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_dir / "data.parquet", index=False, engine="pyarrow")


def test_flatten_one_produces_contract_columns(tmp_path: Path) -> None:
    annotated_dir = tmp_path / "annotated"
    genotypes_dir = tmp_path / "genotypes"
    gene_list = tmp_path / "GeneList.txt"

    gene_list.write_text("symbol\tglobalScore\nTP53\t0.95\n", encoding="utf-8")
    genotypes_dir.mkdir(parents=True, exist_ok=True)
    (genotypes_dir / "samples.txt").write_text("S001\n", encoding="utf-8")

    ann_df = pd.DataFrame(
        [
            {
                "chromosome": "1",
                "position": 101,
                "ref": "A",
                "alt": "G",
                "gencode": {
                    "genes": ["TP53"],
                    "consequence": ["missense_variant"],
                    "region_type": "exonic",
                    "transcripts": ["TP53:NM_000546:exon5:c.215C>G:p.Pro72Arg"],
                },
                "clnsig": "Pathogenic",
                "apc_protein_function_v3": 0.77,
                "bravo_af": 0.02,
            }
        ]
    )
    _write_fixture_parquet(annotated_dir, "1", ann_df)

    geno_df = pd.DataFrame(
        [
            {
                "chromosome": "1",
                "position": 101,
                "ref": "A",
                "alt": "G",
                "dosages": [1.0],
            }
        ]
    )
    _write_fixture_parquet(genotypes_dir, "1", geno_df)

    flat_df, summary = flatten_one(
        annotated_dir=annotated_dir,
        genotypes_dir=genotypes_dir,
        gene_list=gene_list,
        sample_id="S001",
        chromosome="1",
        dosage_threshold=0.0,
    )

    assert len(flat_df) == 1
    assert summary["final_carried_gene_matched_rows"] == 1
    assert flat_df.iloc[0]["matched_gene"] == "TP53"
    assert contract.validate_flat_df(flat_df, raise_on_missing=False) == []
    assert list(flat_df.columns[: len(contract.REQUIRED_FLAT_COLUMNS)]) == contract.REQUIRED_FLAT_COLUMNS
    assert "impact_fallback_flags" in flat_df.columns
    assert "impact_fallback_count" in flat_df.columns
    assert int(flat_df.iloc[0]["impact_fallback_count"]) == 0
    assert int(summary["rows_with_compatibility_fallbacks"]) == 0
    assert summary["compatibility_fallback_counts"] == {}


def test_validate_flat_df_raises_on_missing_required_column() -> None:
    df = pd.DataFrame([
        {
            "sample_id": "S001",
            "variant_id": 1,
            "chromosome": "1",
            "position": 101,
            "ref": "A",
            "alt": "G",
            "allele": "G",
            "dosage": 1.0,
            # maf intentionally missing
            "vid": "1:101:A:G",
            "variant_vcf": "1:101:A:G",
            "matched_gene": "TP53",
            "matched_gene_score": 0.95,
            "genecode_comprehensive_info": "TP53 transcript",
            "genecode_comprehensive_exonic_category": "nonsynonymous SNV",
            "refseq_exonic_category": "nonsynonymous SNV",
            "ucsc_exonic_category": "nonsynonymous SNV",
            "clnsig": "Pathogenic",
            "apc_protein_function_v3": 0.77,
        }
    ])
    with pytest.raises(ValueError, match="missing required columns"):
        contract.validate_flat_df(df)


def test_flatten_one_tracks_fallback_provenance(tmp_path: Path) -> None:
    annotated_dir = tmp_path / "annotated"
    genotypes_dir = tmp_path / "genotypes"
    gene_list = tmp_path / "GeneList.txt"

    gene_list.write_text("symbol\tglobalScore\nTP53\t0.95\n", encoding="utf-8")
    genotypes_dir.mkdir(parents=True, exist_ok=True)
    (genotypes_dir / "samples.txt").write_text("S001\n", encoding="utf-8")

    ann_df = pd.DataFrame(
        [
            {
                "chromosome": "1",
                "position": 101,
                "ref": "A",
                "alt": "G",
                "gencode": {
                    "genes": ["TP53"],
                    "consequence": ["missense_variant"],
                    "region_type": "exonic",
                    # transcripts intentionally missing to trigger fallback
                },
                # clnsig intentionally missing
                # apc intentionally missing
                # af sources intentionally missing
            }
        ]
    )
    _write_fixture_parquet(annotated_dir, "1", ann_df)

    geno_df = pd.DataFrame(
        [
            {
                "chromosome": "1",
                "position": 101,
                "ref": "A",
                "alt": "G",
                "dosages": [1.0],
            }
        ]
    )
    _write_fixture_parquet(genotypes_dir, "1", geno_df)

    flat_df, summary = flatten_one(
        annotated_dir=annotated_dir,
        genotypes_dir=genotypes_dir,
        gene_list=gene_list,
        sample_id="S001",
        chromosome="1",
        dosage_threshold=0.0,
    )

    assert len(flat_df) == 1
    flags = str(flat_df.iloc[0]["impact_fallback_flags"])
    assert "GENECODE_INFO_FALLBACK_TO_MATCHED_GENES" in flags
    assert "CLNSIG_MISSING_DEFAULT_EMPTY" in flags
    assert "APC_PROTEIN_FUNCTION_DEFAULT_0" in flags
    assert "MAF_DEFAULT_0" in flags
    assert int(flat_df.iloc[0]["impact_fallback_count"]) == 4
    assert int(summary["rows_with_compatibility_fallbacks"]) == 1
    assert summary["compatibility_fallback_counts"] == {
        "APC_PROTEIN_FUNCTION_DEFAULT_0": 1,
        "CLNSIG_MISSING_DEFAULT_EMPTY": 1,
        "GENECODE_INFO_FALLBACK_TO_MATCHED_GENES": 1,
        "MAF_DEFAULT_0": 1,
    }


@pytest.mark.parametrize(
    ("consequence", "ref", "alt", "expected"),
    [
        (["inframe_insertion"], "T", "TGGC", "nonframeshift insertion"),
        (["inframe_deletion"], "TGGC", "T", "nonframeshift deletion"),
        (["frameshift_variant"], "T", "TG", "frameshift insertion"),
        (["nonsynonymous_variant"], "A", "G", "nonsynonymous SNV"),
        (["synonymous_variant"], "A", "G", "synonymous SNV"),
        (["stop_lost"], "A", "G", "stoploss"),
    ],
)
def test_map_legacy_exonic_category_avoids_substring_mismatches(
    consequence: list[str], ref: str, alt: str, expected: str
) -> None:
    assert map_legacy_exonic_category(consequence, ref=ref, alt=alt) == expected


@pytest.mark.parametrize("chrom", ["X", "Y"])
def test_flatten_one_supports_sex_chromosomes(tmp_path: Path, chrom: str) -> None:
    annotated_dir = tmp_path / "annotated"
    genotypes_dir = tmp_path / "genotypes"
    gene_list = tmp_path / "GeneList.txt"

    gene_list.write_text("symbol\tglobalScore\nTP53\t0.95\n", encoding="utf-8")
    genotypes_dir.mkdir(parents=True, exist_ok=True)
    (genotypes_dir / "samples.txt").write_text("S001\n", encoding="utf-8")

    ann_df = pd.DataFrame(
        [
            {
                "chromosome": chrom,
                "position": 101,
                "ref": "A",
                "alt": "G",
                "gencode": {
                    "genes": ["TP53"],
                    "consequence": ["missense_variant"],
                    "region_type": "exonic",
                    "transcripts": ["TP53:NM_000546:exon5:c.215C>G:p.Pro72Arg"],
                },
                "clnsig": "Pathogenic",
                "apc_protein_function_v3": 0.77,
                "bravo_af": 0.02,
            }
        ]
    )
    _write_fixture_parquet(annotated_dir, chrom, ann_df)

    geno_df = pd.DataFrame(
        [
            {
                "chromosome": chrom,
                "position": 101,
                "ref": "A",
                "alt": "G",
                "dosages": [1.0],
            }
        ]
    )
    _write_fixture_parquet(genotypes_dir, chrom, geno_df)

    flat_df, summary = flatten_one(
        annotated_dir=annotated_dir,
        genotypes_dir=genotypes_dir,
        gene_list=gene_list,
        sample_id="S001",
        chromosome=chrom,
        dosage_threshold=0.0,
    )

    assert len(flat_df) == 1
    assert str(flat_df.iloc[0]["chromosome"]) == chrom
    assert summary["chromosome"] == chrom


def test_flatten_one_skips_missing_dosage_for_target_sample(tmp_path: Path) -> None:
    annotated_dir = tmp_path / "annotated"
    genotypes_dir = tmp_path / "genotypes"
    gene_list = tmp_path / "GeneList.txt"

    gene_list.write_text("symbol\tglobalScore\nTP53\t0.95\n", encoding="utf-8")
    genotypes_dir.mkdir(parents=True, exist_ok=True)
    # Two samples so the target sample index points to the second element.
    (genotypes_dir / "samples.txt").write_text("S000\nS001\n", encoding="utf-8")

    ann_df = pd.DataFrame(
        [
            {
                "chromosome": "1",
                "position": 101,
                "ref": "A",
                "alt": "G",
                "gencode": {
                    "genes": ["TP53"],
                    "consequence": ["missense_variant"],
                    "region_type": "exonic",
                    "transcripts": ["TP53:NM_000546:exon5:c.215C>G:p.Pro72Arg"],
                },
                "clnsig": "Pathogenic",
                "apc_protein_function_v3": 0.77,
                "bravo_af": 0.02,
            }
        ]
    )
    _write_fixture_parquet(annotated_dir, "1", ann_df)

    geno_df = pd.DataFrame(
        [
            {
                "chromosome": "1",
                "position": 101,
                "ref": "A",
                "alt": "G",
                # second sample dosage is missing; row should be dropped for S001
                "dosages": [1.0, None],
            }
        ]
    )
    _write_fixture_parquet(genotypes_dir, "1", geno_df)

    flat_df, summary = flatten_one(
        annotated_dir=annotated_dir,
        genotypes_dir=genotypes_dir,
        gene_list=gene_list,
        sample_id="S001",
        chromosome="1",
        dosage_threshold=0.0,
    )

    assert len(flat_df) == 0
    assert summary["carried_variants_for_sample"] == 0
    assert summary["final_carried_gene_matched_rows"] == 0


def test_flatten_one_supports_alt_with_multiallelic_style_string(tmp_path: Path) -> None:
    annotated_dir = tmp_path / "annotated"
    genotypes_dir = tmp_path / "genotypes"
    gene_list = tmp_path / "GeneList.txt"

    gene_list.write_text("symbol\tglobalScore\nTP53\t0.95\n", encoding="utf-8")
    genotypes_dir.mkdir(parents=True, exist_ok=True)
    (genotypes_dir / "samples.txt").write_text("S001\n", encoding="utf-8")

    # This simulates a row where ALT carries a comma-separated allele list.
    ann_df = pd.DataFrame(
        [
            {
                "chromosome": "1",
                "position": 202,
                "ref": "C",
                "alt": "A,G",
                "gencode": {
                    "genes": ["TP53"],
                    "consequence": ["missense_variant"],
                    "region_type": "exonic",
                    "transcripts": ["TP53:NM_000546:exon6:c.300C>A"],
                },
                "clnsig": "Likely_pathogenic",
                "apc_protein_function_v3": 0.6,
                "bravo_af": 0.01,
            }
        ]
    )
    _write_fixture_parquet(annotated_dir, "1", ann_df)

    geno_df = pd.DataFrame(
        [
            {
                "chromosome": "1",
                "position": 202,
                "ref": "C",
                "alt": "A,G",
                "dosages": [1.0],
            }
        ]
    )
    _write_fixture_parquet(genotypes_dir, "1", geno_df)

    flat_df, summary = flatten_one(
        annotated_dir=annotated_dir,
        genotypes_dir=genotypes_dir,
        gene_list=gene_list,
        sample_id="S001",
        chromosome="1",
        dosage_threshold=0.0,
    )

    assert len(flat_df) == 1
    assert flat_df.iloc[0]["alt"] == "A,G"
    assert flat_df.iloc[0]["allele"] == "A,G"
    assert summary["final_carried_gene_matched_rows"] == 1
