from __future__ import annotations

import json
from pathlib import Path

from impact_snv.cli import main
from impact_snv.gene_lists import ResolvedPhenotypeTargets, resolve_sample_gene_list, select_search_hit


def test_select_search_hit_prefers_exact_name_match() -> None:
    hit, strategy, rank = select_search_hit(
        "Eczema",
        [
            {"id": "MONDO_1", "name": "Wiskott-Aldrich syndrome"},
            {"id": "HP_0000964", "name": "Eczema"},
        ],
    )
    assert hit["id"] == "HP_0000964"
    assert strategy == "exact_name_search_hit"
    assert rank == 2


def test_resolve_sample_gene_list_uses_manifest_prefix_fallback(tmp_path: Path) -> None:
    gene_list = tmp_path / "Case1.GeneList.txt"
    gene_list.write_text("symbol\tglobalScore\nTP53\t0.95\n", encoding="utf-8")
    manifest = tmp_path / "sample_gene_lists.tsv"
    manifest.write_text(
        "sample_id\tgene_list_path\nCase1\t" + str(gene_list) + "\n",
        encoding="utf-8",
    )
    assert resolve_sample_gene_list("Case1_proband", gene_list_manifest=manifest) == gene_list


def test_build_gene_lists_aggregates_union_and_max_scores(tmp_path: Path, monkeypatch) -> None:
    samples_manifest = tmp_path / "samples.csv"
    samples_manifest.write_text(
        "sample_id,vcf_path,assembly,sex\n"
        "Case1,/tmp/case1.vcf.gz,GRCh38,Female\n",
        encoding="utf-8",
    )
    phenotypes_manifest = tmp_path / "phenotypes.csv"
    phenotypes_manifest.write_text(
        "sample_id,phenotype,hpo_id,status\n"
        "Case1,Eczema,,observed\n"
        "Case1,Atopic dermatitis,,observed\n"
        "Case1,Ignored phenotype,,excluded\n",
        encoding="utf-8",
    )

    class FakeClient:
        def __init__(self, api_url: str):
            self.api_url = api_url

        def resolve_phenotype_targets(
            self,
            *,
            sample_id: str,
            phenotype: str,
            hpo_id: str,
            status: str,
            max_search_hits: int,
            page_size: int,
        ) -> ResolvedPhenotypeTargets:
            if phenotype == "Eczema":
                gene_scores = {"FLG": 0.70, "IL13": 0.55}
                disease_id = "HP_0000964"
            elif phenotype == "Atopic dermatitis":
                gene_scores = {"FLG": 0.91, "STAT6": 0.42}
                disease_id = "EFO_0000274"
            else:  # pragma: no cover - excluded by status filter
                raise AssertionError(phenotype)
            return ResolvedPhenotypeTargets(
                sample_id=sample_id,
                phenotype=phenotype,
                hpo_id=hpo_id,
                status=status,
                query_text=phenotype,
                disease_id=disease_id,
                disease_name=phenotype,
                resolution_strategy="top_search_hit",
                search_hit_rank=1,
                associated_target_count=len(gene_scores),
                gene_count=len(gene_scores),
                max_global_score=max(gene_scores.values()),
                top_gene_symbols=tuple(sorted(gene_scores, key=gene_scores.get, reverse=True)),
                gene_scores=gene_scores,
            )

    monkeypatch.setattr("impact_snv.gene_lists.OpenTargetsClient", FakeClient)

    out_dir = tmp_path / "gene_lists_out"
    rc = main(
        [
            "build-gene-lists",
            "--samples-manifest",
            str(samples_manifest),
            "--phenotypes-manifest",
            str(phenotypes_manifest),
            "--out-dir",
            str(out_dir),
        ]
    )
    assert rc == 0

    gene_list = out_dir / "gene_lists" / "Case1.GeneList.txt"
    assert gene_list.exists()
    assert gene_list.read_text(encoding="utf-8").splitlines() == [
        "symbol\tglobalScore",
        "FLG\t0.91",
        "IL13\t0.55",
        "STAT6\t0.42",
    ]

    manifest_tsv = out_dir / "sample_gene_lists.tsv"
    manifest_json = out_dir / "sample_gene_lists.json"
    assert manifest_tsv.exists()
    assert manifest_json.exists()
    payload = json.loads(manifest_json.read_text(encoding="utf-8"))
    assert payload["command"] == "impact-snv build-gene-lists"
    assert payload["sample_count"] == 1
    assert payload["phenotype_query_count"] == 2
    assert payload["sample_gene_lists"][0]["gene_count"] == 3