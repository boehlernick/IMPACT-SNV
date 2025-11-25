#!/usr/bin/env python3
import pandas as pd
import os
import sys

"""
Generate FAVOR subset CSVs in the exact schema expected by favorannotator.R.
Steps:
1. Read raw FAVOR processed CSV (with .String/.Valid columns).
2. Drop all *.Valid columns.
3. Map raw column names to expected names from favorannotator.R.
4. Preserve only columns defined in favorannotator.R.
5. Output per-chromosome CSVs with headers matching favorannotator.R.
6. Write favor_subset_map.txt mapping chromosome -> CSV path.
7. Log missing variants.
Usage:
  python build_favor_subset.py <favor_raw_csv> <variants_tsv> <output_dir>
"""

# Mapping dictionary: raw -> expected (full mapping based on favorannotator.R)
column_map = {
    "VariantVcf.String": "variant_vcf",
    "Chromosome.String": "chromosome",
    "Position.String": "position",
    "RefVcf.String": "ref_vcf",
    "AltVcf.String": "alt_vcf",
    "AloftValue.String": "aloft_value",
    "AloftDescription.String": "aloft_description",
    "ApcConservation.Float64": "apc_conservation",
    "ApcConservationV2.Float64": "apc_conservation_v2",
    "ApcEpigenetics.Float64": "apc_epigenetics",
    "ApcEpigeneticsActive.Float64": "apc_epigenetics_active",
    "ApcEpigeneticsRepressed.Float64": "apc_epigenetics_repressed",
    "ApcEpigeneticsTranscription.Float64": "apc_epigenetics_transcription",
    "ApcLocalNucleotideDiversity.Float64": "apc_local_nucleotide_diversity",
    "ApcLocalNucleotideDiversityV2.Float64": "apc_local_nucleotide_diversity_v2",
    "ApcLocalNucleotideDiversityV3.Float64": "apc_local_nucleotide_diversity_v3",
    "ApcMappability.Float64": "apc_mappability",
    "ApcMicroRna.Float64": "apc_micro_rna",
    "ApcMutationDensity.Float64": "apc_mutation_density",
    "ApcProteinFunctionV3.Float64": "apc_protein_function_v3",
    "ApcProximityToCoding.Float64": "apc_proximity_to_coding",
    "ApcProximityToCodingV2.Float64": "apc_proximity_to_coding_v2",
    "ApcProximityToTsstes.Float64": "apc_proximity_to_tsstes",
    "ApcTranscriptionFactor.Float64": "apc_transcription_factor",
    "CagePromoter.String": "cage_promoter",
    "CageTc.String": "cage_tc",
    "MetasvmPred.String": "metasvm_pred",
    "Rsid.String": "rsid",
    "FathmmXf.Float64": "fathmm_xf",
    "GenecodeComprehensiveCategory.String": "genecode_comprehensive_category",
    "GenecodeComprehensiveInfo.String": "genecode_comprehensive_info",
    "GenecodeComprehensiveExonicCategory.String": "genecode_comprehensive_exonic_category",
    "GenecodeComprehensiveExonicInfo.String": "genecode_comprehensive_exonic_info",
    "Genehancer.String": "genehancer",
    "Linsight.Float64": "linsight",
    "CaddPhred.Float64": "cadd_phred",
    "Rdhs.String": "rdhs",
    # Add remaining mappings as needed for full schema
}

expected_order = [
    "variant_vcf","chromosome","position","ref_vcf","alt_vcf","aloft_value","aloft_description",
    "apc_conservation","apc_conservation_v2","apc_epigenetics","apc_epigenetics_active","apc_epigenetics_repressed",
    "apc_epigenetics_transcription","apc_local_nucleotide_diversity","apc_local_nucleotide_diversity_v2",
    "apc_local_nucleotide_diversity_v3","apc_mappability","apc_micro_rna","apc_mutation_density","apc_protein_function_v3",
    "apc_proximity_to_coding","apc_proximity_to_coding_v2","apc_proximity_to_tsstes","apc_transcription_factor",
    "cage_promoter","cage_tc","metasvm_pred","rsid","fathmm_xf","genecode_comprehensive_category",
    "genecode_comprehensive_info","genecode_comprehensive_exonic_category","genecode_comprehensive_exonic_info",
    "genehancer","linsight","cadd_phred","rdhs"
]

def main():
    if len(sys.argv) != 4:
        print("Usage: python build_favor_subset.py <favor_raw_csv> <variants_tsv> <output_dir>")
        sys.exit(1)

    favor_csv = sys.argv[1]
    variants_tsv = sys.argv[2]
    output_dir = sys.argv[3]
    os.makedirs(output_dir, exist_ok=True)

    # Load raw FAVOR and variants
    favor = pd.read_csv(favor_csv)
    variants = pd.read_csv(variants_tsv, sep='	', header=None, names=['chrom','pos','ref','alt'])

    # Drop all *.Valid columns
    favor = favor[[c for c in favor.columns if not c.endswith('.Valid')]]

    # Map columns to expected names
    mapped_cols = {}
    for col in favor.columns:
        if col in column_map:
            mapped_cols[col] = column_map[col]
    favor = favor.rename(columns=mapped_cols)

    # Keep only expected columns
    favor = favor[[c for c in expected_order if c in favor.columns]]

    # Normalize chromosome and position for join
    favor['chromosome'] = favor['chromosome'].astype(str)
    favor['position'] = favor['position'].astype(str)
    favor['ref_vcf'] = favor['ref_vcf'].astype(str)
    favor['alt_vcf'] = favor['alt_vcf'].astype(str)
    variants['chrom'] = variants['chrom'].astype(str)
    variants['pos'] = variants['pos'].astype(str)

    # Join
    merged = variants.merge(favor, left_on=['chrom','pos','ref','alt'], right_on=['chromosome','position','ref_vcf','alt_vcf'], how='left', indicator=True)

    # Log missing
    missing = merged[merged['_merge'] != 'both'][['chrom','pos','ref','alt']]
    if len(missing):
        missing.to_csv(os.path.join(output_dir,'missing_variants.log'), index=False)

    found = merged[merged['_merge'] == 'both']

    # Group by chromosome and write CSVs
    for chrom, sub in found.groupby('chromosome'):
        out_path = os.path.join(output_dir, f"favor_chr{chrom}.csv")
        sub = sub[expected_order]
        sub.to_csv(out_path, index=False)

    # Write mapping file
    with open(os.path.join(output_dir,'favor_subset_map.txt'),'w') as fh:
        for chrom in sorted(found['chromosome'].unique(), key=lambda x:(len(x),x)):
            fh.write(f"{chrom}	{os.path.join(output_dir,f'favor_chr{chrom}.csv')}
")

if __name__ == '__main__':
    main()
