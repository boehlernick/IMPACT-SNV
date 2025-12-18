<!-- dx-header -->
# IMPACT SNV and InDel vcf processing (DNAnexus Platform App)

Merges SNV and InDel vcf input files into a single merged_vcf.gz output file

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.
<!-- /dx-header -->

## Overview

The IMPACT SNV and InDel VCF Processing app is designed to merge multiple VCF or VCF.GZ files containing Single Nucleotide Variants (SNVs) and Insertions/Deletions (InDels) into a single, compressed VCF.GZ file. This app leverages BCFtools for efficient merging and formatting of the input files, ensuring compatibility with downstream analysis pipelines.

## Features

- **Input**: Accepts multiple VCF or VCF.GZ files.
- **Output**: Produces a single merged VCF.GZ file.
- **Automated Installation**: Automatically installs BCFtools if not already present.
- **Chromosome Separation**: Merges and formats VCF files by chromosome.

## Usage

### Inputs

- **input_vcfs**: An array of VCF or VCF.GZ files to be merged.

### Outputs

- **split_merged_vcf**: merged sample VCF.gz files for chromosomes 1-22