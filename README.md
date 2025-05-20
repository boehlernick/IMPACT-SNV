# IMPACT-SNV
This repository contains the full IMPACT-SNV branch of the IMPACT pipeline, which processes and prioritizes single nucleotide variants (SNVs) and indels for rare disease analysis using the FAVOR database. 

## Pipeline Overview

The pipeline consists of four main steps, each available as a DNAnexus applet:

1. **Step1_vcf2gds**: Converts VCF files to GDS format for efficient downstream processing.
2. **Step2_vcf_merge**: Merges multiple VCF or VCF.GZ files into a single, chromosome-separated VCF.GZ file.
3. **Step3_favorannotator-rap**: Functionally annotates GDS files using the FAVOR database, producing annotated GDS (AGDS) files.
4. **Step4_impact_prioritization**: Scores and prioritizes variants using gene-disease association data and functional scores. Outputs per-sample SNV_IMPACT.gds files. Now available as a DNAnexus applet with robust scoring, tiering, and merging features.

Each step is contained in its own subdirectory with detailed documentation and scripts.

## Requirements
- R (with Bioconductor packages: SeqArray, SeqVarTools, etc.)
- Python 3 (for merging scripts)
- BCFtools, bgzip, tabix (for VCF processing)

Note that these tools were made to be deployed on the DNAnexus Research Analysis platform.

See each step’s README for specific dependencies.

## Usage

### Step 1: VCF Merge
See `step1_vcf_merge/Readme.md` for details. Example command:
```sh
dx run vcf_merge -ivcfs=input1.vcf.gz,input2.vcf.gz -o merged.vcf.gz
```

### Step 2: VCF to GDS Conversion
See `step2_vcf2gds/README.md` for details. Example command:
```sh
dx run /path/to/install/apps/vcf2gds \
  -ivcf_file=/path/to/vcf/file/to/convert/my.vcf.gz \
  -igds_filename=my.gds \
  --priority high \
  -y
```

### Step 3: Functional Annotation
See `step3_favorannotator-rap/README.md` for details. Example command:
```sh
dx run favorannotator -igds=input.gds -o favor_merged_chr*.gds
```

### Step 4: Variant Prioritization (New DNAnexus Applet)
See `step4_impact_prioritization/README.md` for details. Example command:
```sh
dx run step4_impact_prioritization -igda=GeneList.txt -igds_files=favor_merged_chr1.gds,favor_merged_chr2.gds,...
```
- `--genelist` (optional): Path to gene-disease association file (default: `GeneList.txt`)
- `--outprefix` (optional): Prefix for output GDS files (default: `anno_merged_`)
- Outputs: All `<sample_id>_SNV_IMPACT.gds` files

## References
- Original tools and documentation:
  - [vcf2gds](https://github.com/drarwood/vcf2gds)
  - [favorannotator-rap](https://github.com/li-lab-genetics/favorannotator-rap)
- For more information, see the subdirectory READMEs.

## Citation
If you use this pipeline, please cite the relevant tools and the IMPACT pipeline publication [pending]

