#!/usr/bin/env python3
"""Sample-aware Open Targets gene list generation and resolution helpers."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
import re
from typing import Any, Optional, Sequence
import urllib.error
import urllib.request

from impact_snv import __version__ as VERSION


DEFAULT_OPENTARGETS_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"
DEFAULT_MAX_SEARCH_HITS = 10
DEFAULT_PAGE_SIZE = 500
DEFAULT_INCLUDED_STATUSES = frozenset({"", "observed"})


SEARCH_DISEASE_QUERY = """
query SearchDisease($term: String!, $size: Int!) {
  search(queryString: $term, entityNames: ["disease"], page: {index: 0, size: $size}) {
    hits {
      id
      name
      entity
      description
    }
  }
}
"""


ASSOCIATED_TARGETS_QUERY = """
query DiseaseAssociatedTargets($diseaseId: String!, $pageIndex: Int!, $pageSize: Int!) {
  disease(efoId: $diseaseId) {
    id
    name
    associatedTargets(page: {index: $pageIndex, size: $pageSize}) {
      count
      rows {
        score
        target {
          id
          approvedSymbol
        }
      }
    }
  }
}
"""


class OpenTargetsError(RuntimeError):
    """Raised when Open Targets querying or phenotype resolution fails."""


@dataclass(frozen=True)
class SampleManifestRecord:
    sample_id: str
    vcf_path: str
    assembly: str
    sex: str


@dataclass(frozen=True)
class PhenotypeManifestRecord:
    sample_id: str
    phenotype: str
    hpo_id: str
    status: str


@dataclass(frozen=True)
class ResolvedPhenotypeTargets:
    sample_id: str
    phenotype: str
    hpo_id: str
    status: str
    query_text: str
    disease_id: str
    disease_name: str
    resolution_strategy: str
    search_hit_rank: int
    associated_target_count: int
    gene_count: int
    max_global_score: Optional[float]
    top_gene_symbols: tuple[str, ...]
    gene_scores: dict[str, float]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class SampleGeneListRecord:
    sample_id: str
    gene_list_path: str
    phenotype_count: int
    gene_count: int
    max_global_score: Optional[float]
    resolution_strategies: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class BuildGeneListsResult:
    command: str
    version: str
    created_utc: str
    samples_manifest: str
    phenotypes_manifest: str
    out_dir: str
    manifest_tsv: str
    manifest_json: str
    api_url: str
    sample_count: int
    phenotype_query_count: int
    sample_gene_lists: tuple[SampleGeneListRecord, ...]
    phenotype_queries: tuple[ResolvedPhenotypeTargets, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            **asdict(self),
            "sample_gene_lists": [record.to_dict() for record in self.sample_gene_lists],
            "phenotype_queries": [record.to_dict() for record in self.phenotype_queries],
        }


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize_token(value: str, fallback: str = "sample") -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    token = re.sub(r"_+", "_", token).strip("._-")
    return token or fallback


def normalize_match_text(value: str) -> str:
    text = str(value or "").lower().strip()
    text = text.replace("_", " ")
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def normalize_ontology_id(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    return text.replace(":", "_")


def detect_delimiter(path: Path) -> str:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
    if "\t" in sample and "," not in sample.splitlines()[0]:
        return "\t"
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except csv.Error:
        return ","


def read_tabular_rows(path: Path) -> list[dict[str, str]]:
    delimiter = detect_delimiter(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Manifest has no header row: {path}")
        return [{key: (value or "") for key, value in row.items()} for row in reader]


def require_columns(path: Path, rows: Sequence[dict[str, str]], required: set[str]) -> None:
    if not rows:
        raise ValueError(f"Manifest is empty: {path}")
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Manifest {path} is missing required columns: {sorted(missing)}")


def read_samples_manifest(path: Path) -> list[SampleManifestRecord]:
    rows = read_tabular_rows(path)
    require_columns(path, rows, {"sample_id"})
    seen: set[str] = set()
    out: list[SampleManifestRecord] = []
    for row in rows:
        sample_id = str(row.get("sample_id", "")).strip()
        if not sample_id:
            raise ValueError(f"Blank sample_id in samples manifest: {path}")
        if sample_id in seen:
            raise ValueError(f"Duplicate sample_id in samples manifest: {sample_id}")
        seen.add(sample_id)
        out.append(
            SampleManifestRecord(
                sample_id=sample_id,
                vcf_path=str(row.get("vcf_path", "")).strip(),
                assembly=str(row.get("assembly", "")).strip(),
                sex=str(row.get("sex", "")).strip(),
            )
        )
    return out


def read_phenotypes_manifest(
    path: Path,
    *,
    included_statuses: Sequence[str] = tuple(DEFAULT_INCLUDED_STATUSES),
) -> list[PhenotypeManifestRecord]:
    rows = read_tabular_rows(path)
    require_columns(path, rows, {"sample_id", "phenotype"})
    allowed = {normalize_match_text(value) for value in included_statuses}
    out: list[PhenotypeManifestRecord] = []
    for row in rows:
        status = str(row.get("status", "")).strip()
        if normalize_match_text(status) not in allowed:
            continue
        sample_id = str(row.get("sample_id", "")).strip()
        phenotype = str(row.get("phenotype", "")).strip()
        if not sample_id or not phenotype:
            raise ValueError(f"Phenotype manifest row is missing sample_id or phenotype: {row}")
        out.append(
            PhenotypeManifestRecord(
                sample_id=sample_id,
                phenotype=phenotype,
                hpo_id=normalize_ontology_id(str(row.get("hpo_id", "")).strip()),
                status=status,
            )
        )
    if not out:
        raise ValueError(f"No phenotype rows remained after status filtering: {path}")
    return out


def group_phenotypes_by_sample(
    samples: Sequence[SampleManifestRecord],
    phenotypes: Sequence[PhenotypeManifestRecord],
) -> dict[str, list[PhenotypeManifestRecord]]:
    sample_ids = {sample.sample_id for sample in samples}
    out: dict[str, list[PhenotypeManifestRecord]] = {sample.sample_id: [] for sample in samples}
    extras = sorted({record.sample_id for record in phenotypes if record.sample_id not in sample_ids})
    if extras:
        raise ValueError(f"Phenotype manifest contains sample_id values not present in samples manifest: {extras}")
    for record in phenotypes:
        out[record.sample_id].append(record)
    missing = sorted([sample_id for sample_id, rows in out.items() if not rows])
    if missing:
        raise ValueError(f"Samples manifest entries are missing phenotypes: {missing}")
    return out


def numeric_score(value: Any) -> Optional[float]:
    try:
        out = float(value)
    except Exception:
        return None
    return out


def max_or_none(values: Sequence[float]) -> Optional[float]:
    return max(values) if values else None


def select_search_hit(query_text: str, hits: Sequence[dict[str, Any]]) -> tuple[dict[str, Any], str, int]:
    if not hits:
        raise OpenTargetsError(f"Open Targets search returned no disease hits for phenotype query: {query_text!r}")
    normalized_query = normalize_match_text(query_text)
    for index, hit in enumerate(hits, start=1):
        if normalize_match_text(str(hit.get("name", ""))) == normalized_query:
            return hit, "exact_name_search_hit", index
    return hits[0], "top_search_hit", 1


def extract_gene_scores(rows: Sequence[dict[str, Any]]) -> dict[str, float]:
    out: dict[str, float] = {}
    for row in rows:
        target = row.get("target") or {}
        symbol = str(target.get("approvedSymbol", "")).strip()
        score = numeric_score(row.get("score"))
        if not symbol or score is None:
            continue
        current = out.get(symbol)
        if current is None or score > current:
            out[symbol] = score
    return out


def sort_gene_scores(gene_scores: dict[str, float]) -> list[tuple[str, float]]:
    return sorted(gene_scores.items(), key=lambda item: (-item[1], item[0]))


class OpenTargetsClient:
    def __init__(self, api_url: str = DEFAULT_OPENTARGETS_API_URL):
        self.api_url = api_url

    def graphql(self, query: str, variables: dict[str, Any]) -> dict[str, Any]:
        payload = json.dumps({"query": query, "variables": variables}).encode("utf-8")
        request = urllib.request.Request(
            self.api_url,
            data=payload,
            headers={"Content-Type": "application/json"},
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                data = json.load(response)
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            raise OpenTargetsError(
                f"Open Targets API HTTP {exc.code} for {self.api_url}: {detail[:2000]}"
            ) from exc
        except urllib.error.URLError as exc:
            raise OpenTargetsError(f"Open Targets API request failed: {exc}") from exc
        errors = data.get("errors") or []
        if errors:
            raise OpenTargetsError(json.dumps(errors, indent=2, sort_keys=True))
        return data.get("data") or {}

    def search_disease(self, query_text: str, *, max_hits: int) -> list[dict[str, Any]]:
        data = self.graphql(SEARCH_DISEASE_QUERY, {"term": query_text, "size": max_hits})
        return list(((data.get("search") or {}).get("hits") or []))

    def associated_targets_page(
        self,
        disease_id: str,
        *,
        page_index: int,
        page_size: int,
    ) -> tuple[str, str, int, list[dict[str, Any]]]:
        data = self.graphql(
            ASSOCIATED_TARGETS_QUERY,
            {"diseaseId": disease_id, "pageIndex": page_index, "pageSize": page_size},
        )
        disease = data.get("disease")
        if disease is None:
            raise OpenTargetsError(f"Open Targets returned no disease entity for ID: {disease_id}")
        associated = disease.get("associatedTargets") or {}
        return (
            str(disease.get("id") or disease_id),
            str(disease.get("name") or ""),
            int(associated.get("count") or 0),
            list(associated.get("rows") or []),
        )

    def fetch_associated_targets(self, disease_id: str, *, page_size: int) -> tuple[str, dict[str, float], int]:
        total_count = 0
        disease_name = ""
        rows: list[dict[str, Any]] = []
        page_index = 0
        while True:
            _, disease_name, page_count, page_rows = self.associated_targets_page(
                disease_id,
                page_index=page_index,
                page_size=page_size,
            )
            total_count = page_count
            rows.extend(page_rows)
            if not page_rows or len(rows) >= page_count or len(page_rows) < page_size:
                break
            page_index += 1
        return disease_name, extract_gene_scores(rows), total_count

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
        disease_id = normalize_ontology_id(hpo_id)
        disease_name = ""
        resolution_strategy = "provided_ontology_id"
        search_hit_rank = 0
        query_text = disease_id or phenotype

        if not disease_id:
            hit, resolution_strategy, search_hit_rank = select_search_hit(
                phenotype,
                self.search_disease(phenotype, max_hits=max_search_hits),
            )
            disease_id = str(hit.get("id") or "").strip()
            disease_name = str(hit.get("name") or "").strip()
            query_text = phenotype
        if not disease_id:
            raise OpenTargetsError(f"Could not resolve phenotype to an Open Targets disease ID: {phenotype!r}")

        fetched_name, gene_scores, associated_target_count = self.fetch_associated_targets(
            disease_id,
            page_size=page_size,
        )
        if not disease_name:
            disease_name = fetched_name
        ordered_genes = sort_gene_scores(gene_scores)
        max_global_score = ordered_genes[0][1] if ordered_genes else None
        return ResolvedPhenotypeTargets(
            sample_id=sample_id,
            phenotype=phenotype,
            hpo_id=normalize_ontology_id(hpo_id),
            status=status,
            query_text=query_text,
            disease_id=disease_id,
            disease_name=disease_name,
            resolution_strategy=resolution_strategy,
            search_hit_rank=search_hit_rank,
            associated_target_count=associated_target_count,
            gene_count=len(gene_scores),
            max_global_score=max_global_score,
            top_gene_symbols=tuple(symbol for symbol, _ in ordered_genes[:10]),
            gene_scores=gene_scores,
        )


def write_gene_list(path: Path, gene_scores: dict[str, float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["symbol\tglobalScore"]
    for symbol, score in sort_gene_scores(gene_scores):
        lines.append(f"{symbol}\t{score:.12g}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest_tsv(path: Path, records: Sequence[SampleGeneListRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["sample_id\tgene_list_path\tphenotype_count\tgene_count\tmax_global_score\tresolution_strategies"]
    for record in records:
        lines.append(
            "\t".join(
                [
                    record.sample_id,
                    record.gene_list_path,
                    str(record.phenotype_count),
                    str(record.gene_count),
                    "" if record.max_global_score is None else f"{record.max_global_score:.12g}",
                    ";".join(record.resolution_strategies),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest_json(path: Path, result: BuildGeneListsResult) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result.to_dict(), indent=2, sort_keys=True), encoding="utf-8")


def build_sample_gene_lists(
    *,
    samples_manifest: Path,
    phenotypes_manifest: Path,
    out_dir: Path,
    client: OpenTargetsClient,
    force: bool = False,
    manifest_tsv: Optional[Path] = None,
    manifest_json: Optional[Path] = None,
    max_search_hits: int = DEFAULT_MAX_SEARCH_HITS,
    page_size: int = DEFAULT_PAGE_SIZE,
    included_statuses: Sequence[str] = tuple(DEFAULT_INCLUDED_STATUSES),
) -> BuildGeneListsResult:
    samples = read_samples_manifest(samples_manifest)
    phenotypes = read_phenotypes_manifest(phenotypes_manifest, included_statuses=included_statuses)
    phenotypes_by_sample = group_phenotypes_by_sample(samples, phenotypes)

    out_dir.mkdir(parents=True, exist_ok=True)
    gene_lists_dir = out_dir / "gene_lists"
    tsv_path = (manifest_tsv or (out_dir / "sample_gene_lists.tsv")).resolve()
    json_path = (manifest_json or (out_dir / "sample_gene_lists.json")).resolve()

    sample_results: list[SampleGeneListRecord] = []
    phenotype_queries: list[ResolvedPhenotypeTargets] = []

    for sample in samples:
        sample_gene_scores: dict[str, float] = {}
        sample_queries: list[ResolvedPhenotypeTargets] = []
        for phenotype_record in phenotypes_by_sample[sample.sample_id]:
            resolved = client.resolve_phenotype_targets(
                sample_id=sample.sample_id,
                phenotype=phenotype_record.phenotype,
                hpo_id=phenotype_record.hpo_id,
                status=phenotype_record.status,
                max_search_hits=max_search_hits,
                page_size=page_size,
            )
            sample_queries.append(resolved)
            phenotype_queries.append(resolved)
            for symbol, score in resolved.gene_scores.items():
                current = sample_gene_scores.get(symbol)
                if current is None or score > current:
                    sample_gene_scores[symbol] = score

        if not sample_gene_scores:
            raise OpenTargetsError(
                f"Open Targets queries produced no gene associations for sample {sample.sample_id!r}"
            )

        gene_list_path = (gene_lists_dir / f"{normalize_token(sample.sample_id)}.GeneList.txt").resolve()
        if gene_list_path.exists() and not force:
            raise FileExistsError(f"Refusing to overwrite existing gene list: {gene_list_path}")
        write_gene_list(gene_list_path, sample_gene_scores)

        sample_results.append(
            SampleGeneListRecord(
                sample_id=sample.sample_id,
                gene_list_path=str(gene_list_path),
                phenotype_count=len(sample_queries),
                gene_count=len(sample_gene_scores),
                max_global_score=max_or_none(list(sample_gene_scores.values())),
                resolution_strategies=tuple(sorted({query.resolution_strategy for query in sample_queries})),
            )
        )

    result = BuildGeneListsResult(
        command="impact-snv build-gene-lists",
        version=VERSION,
        created_utc=utc_now_iso(),
        samples_manifest=str(samples_manifest.resolve()),
        phenotypes_manifest=str(phenotypes_manifest.resolve()),
        out_dir=str(out_dir.resolve()),
        manifest_tsv=str(tsv_path),
        manifest_json=str(json_path),
        api_url=client.api_url,
        sample_count=len(sample_results),
        phenotype_query_count=len(phenotype_queries),
        sample_gene_lists=tuple(sample_results),
        phenotype_queries=tuple(phenotype_queries),
    )
    write_manifest_tsv(tsv_path, sample_results)
    write_manifest_json(json_path, result)
    return result


def manifest_record_path(path_text: str, manifest_path: Path) -> Path:
    path = Path(path_text)
    if not path.is_absolute():
        path = (manifest_path.parent / path).resolve()
    return path


@lru_cache(maxsize=16)
def load_gene_list_manifest(path: Path | str) -> dict[str, Path]:
    manifest_path = Path(path).resolve()
    if manifest_path.suffix.lower() == ".json":
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        records = payload.get("sample_gene_lists") or payload.get("samples") or []
    else:
        records = read_tabular_rows(manifest_path)
    mapping: dict[str, Path] = {}
    for record in records:
        sample_id = str(record.get("sample_id", "")).strip()
        gene_list_path = str(record.get("gene_list_path") or record.get("gene_list") or "").strip()
        if not sample_id or not gene_list_path:
            continue
        resolved_path = manifest_record_path(gene_list_path, manifest_path)
        if not resolved_path.exists():
            raise FileNotFoundError(
                f"Gene list manifest entry for {sample_id!r} points to a missing file: {resolved_path}"
            )
        mapping[sample_id] = resolved_path
    if not mapping:
        raise ValueError(f"Gene list manifest contained no usable sample mappings: {manifest_path}")
    return mapping


def resolve_sample_gene_list(
    sample_id: str,
    *,
    gene_list: Optional[Path] = None,
    gene_list_manifest: Optional[Path] = None,
) -> Path:
    if gene_list_manifest is not None:
        mapping = load_gene_list_manifest(gene_list_manifest)
        exact = mapping.get(sample_id)
        if exact is not None:
            return exact
        prefix_matches = [
            (key, path)
            for key, path in mapping.items()
            if sample_id.startswith(f"{key}_")
        ]
        if prefix_matches:
            prefix_matches.sort(key=lambda item: len(item[0]), reverse=True)
            return prefix_matches[0][1]
        if gene_list is None:
            known = ", ".join(sorted(mapping.keys())[:20])
            raise KeyError(
                f"No gene list mapping found for sample {sample_id!r} in {gene_list_manifest}. "
                f"Known manifest sample_id values include: {known}"
            )
    if gene_list is None:
        raise ValueError("A gene list path or gene list manifest is required")
    resolved = Path(gene_list).resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"Gene list not found: {resolved}")
    return resolved


def run_build_gene_lists(args: Any) -> int:
    out_dir = Path(args.out_dir).expanduser().resolve()
    samples_manifest = Path(args.samples_manifest).expanduser().resolve()
    phenotypes_manifest = Path(args.phenotypes_manifest).expanduser().resolve()
    manifest_tsv = Path(args.manifest_tsv).expanduser().resolve() if getattr(args, "manifest_tsv", None) else None
    manifest_json = Path(args.manifest_json).expanduser().resolve() if getattr(args, "manifest_json", None) else None
    client = OpenTargetsClient(api_url=str(getattr(args, "api_url", DEFAULT_OPENTARGETS_API_URL)))
    result = build_sample_gene_lists(
        samples_manifest=samples_manifest,
        phenotypes_manifest=phenotypes_manifest,
        out_dir=out_dir,
        client=client,
        force=bool(getattr(args, "force", False)),
        manifest_tsv=manifest_tsv,
        manifest_json=manifest_json,
        max_search_hits=int(getattr(args, "max_search_hits", DEFAULT_MAX_SEARCH_HITS) or DEFAULT_MAX_SEARCH_HITS),
        page_size=int(getattr(args, "page_size", DEFAULT_PAGE_SIZE) or DEFAULT_PAGE_SIZE),
    )
    if not bool(getattr(args, "quiet", False)):
        print("Sample-specific GeneList generation completed successfully")
        print(f"  samples: {result.sample_count}")
        print(f"  phenotype_queries: {result.phenotype_query_count}")
        print(f"  out_dir: {result.out_dir}")
        print(f"  manifest_tsv: {result.manifest_tsv}")
        print(f"  manifest_json: {result.manifest_json}")
        for record in result.sample_gene_lists:
            print(
                f"  {record.sample_id}: {record.gene_count} genes from {record.phenotype_count} phenotype(s)"
            )
    return 0


__all__ = [
    "BuildGeneListsResult",
    "DEFAULT_OPENTARGETS_API_URL",
    "OpenTargetsClient",
    "OpenTargetsError",
    "ResolvedPhenotypeTargets",
    "SampleGeneListRecord",
    "build_sample_gene_lists",
    "load_gene_list_manifest",
    "resolve_sample_gene_list",
    "run_build_gene_lists",
    "select_search_hit",
]