# Synthetic GRCh38 VCF Test Suite

This folder contains three small synthetic VCF files for testing variant annotation and prioritization pipelines.

## Files

- `synthetic_skeletal_female.vcf` — sample `SYNTH_SKEL_F`; chromosomes present: 1-8 and X; high-priority truth variant: FGFR3 c.1138G>A / p.Gly380Arg.
- `synthetic_cf_male.vcf` — sample `SYNTH_CF_M`; chromosomes present: 7, 9-16, and Y; high-priority truth variant: CFTR F508del.
- `synthetic_cancer_mito_female.vcf` — sample `SYNTH_CA_MT_F`; chromosomes present: 17-22, X, and MT; high-priority truth variants: BRCA1 c.181T>G and MT-TL1 m.3243A>G.
- `synthetic_vcf_manifest.tsv` — machine-readable expected annotations and prioritization outcomes.
- `validation_commands.sh` — strict validation commands; intentionally does not use `bcftools norm -c s`.

## Design choices

- Chromosome names use `1-22`, `X`, `Y`, and `MT` without `chr` prefixes.
- Headers include GRCh38 primary chromosome lengths.
- Non-clinical coverage controls are synthetic, use `ID=.` and `GENE=SYNTHETIC_INTERGENIC`, and are placed at position 1000 with `REF=N` to target the start-of-contig gap/N region common to GRCh38 reference chromosomes. Validate these against your exact GRCh38 FASTA with `bcftools norm -f GRCh38.fa -c e`.
- The mitochondrial pathogenic variant includes `FORMAT/AF` to represent heteroplasmy-like allele fraction while retaining `GT` for compatibility.

## Expected edge cases

- Chromosomes absent from individual files while covered across the suite.
- Haploid Y genotype: `GT=1` with biallelic `AD=0,30`.
- Mitochondrial heteroplasmy proxy: `AD=80,40; AF=0.333`.
- Multiallelic parser test: chr19 `N>A,C`.
- Small insertion parser test: chr22 `N>NA`.
- Low-quality filter control: chr9 `FILTER=LowQual`.
