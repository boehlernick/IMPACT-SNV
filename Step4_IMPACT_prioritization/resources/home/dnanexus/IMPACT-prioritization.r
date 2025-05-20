
#!/usr/bin/env Rscript

# Set CRAN mirror for non-interactive environments
options(repos = c(CRAN = "https://cloud.r-project.org"))

# IMPACT-SNV Variant Prioritization Script
# Processes annotated GDS files and assigns pathogenicity scores and
# tiers to variants.
# Usage: Rscript IMPACT-prioritization.r --genelist GeneList.txt
#   --outprefix anno_merged_

# Argument parsing
suppressPackageStartupMessages({
if (!requireNamespace("optparse", quietly = TRUE))
install.packages("optparse")
library(optparse)
})

option_list <- list(
make_option(c("-g", "--genelist"), type = "character",
  default = "GeneList.txt",
  help = "Gene-disease association file [default %default]"),
make_option(c("-o", "--outprefix"), type = "character",
  default = "anno_merged_",
  help = "Output GDS file prefix [default %default]"),
make_option(c("-p", "--prefix"), type = "character",
  default = "merged_chr",
  help = "Prefix for input GDS files [default %default]"),
make_option(c("--pattern"), type = "character",
  default = "merged_chr(.*)\\.gds",
  help = "Regex pattern to extract chromosome from filenames [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
gda_file <- opt$genelist
outprefix <- opt$outprefix
prefix <- opt$prefix
pattern <- opt$pattern

# Check for matching files before loading packages

gds_files <- list.files(
  pattern = "^merged_chr.*\\.gds$",
  full.names = TRUE
)
if (length(gds_files) == 0) {
stop(paste("No GDS files found matching prefix:", prefix))
}

genelist_path <- list.files(
  pattern = "GeneList.txt$",
  full.names = TRUE
)

if (length(genelist_path) == 0) {
stop("GeneList.txt not found in input directory.")
}
Open_Target_data <- read.table(genelist_path[1], sep = "\t", header = TRUE)

# Load packages only if files are found
suppressPackageStartupMessages({
if (!requireNamespace("rlang", quietly = TRUE)) install.packages("rlang")
if (!requireNamespace("cli", quietly = TRUE)) install.packages("cli")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SeqArray", "SeqVarTools"))
library(rlang)
library(cli)
library(dplyr)
library(stringr)
library(parallel)
library(readr)
library(SeqArray)
library(SeqVarTools)
library(tidyr)
})

# Utility functions
extract_symbols <- function(entry) {
if (is.na(entry) || entry %in% c("", "NONE", "NONE(dist=NONE)")) return(NULL)
entry <- gsub("\\([^d][^i][^s][^t][^=]*\\)", "", entry)
symbols <- unlist(strsplit(entry, ","))
symbols <- gsub("\\(dist=.*\\)", "", symbols)
trimws(symbols)
}
matches_criteria <- function(symbol, cleaned_entry) {
symbol <- gsub("([()])", "\\1", symbol)
pattern <- paste0(symbol, "(,|\\(dist=[0-9]{1,4}\\)|$)")
grepl(pattern, cleaned_entry)
}
get_exonic_indices <- function(aGDS, category, values) {
unlist(lapply(values, function(x) which(seqGetData(aGDS, category) == x)))
}

score_variants <- function(aGDS, Open_Target_data, outprefix, chr, gdsfile) {
  all_variants <- seqGetData(aGDS, "variant.id")
  clnsig <- seqGetData(aGDS, "annotation/info/FunctionalAnnotation/clnsig")
  genecode_info <- seqGetData(aGDS, "annotation/info/FunctionalAnnotation/genecode_comprehensive_info")
  global_score_vector <- numeric(length(genecode_info))
  matched_symbols <- character()
  for (i in seq_along(genecode_info)) {
    entry <- genecode_info[i]
    symbols <- extract_symbols(entry)
    if (!is.null(symbols)) {
      cleaned_entry <- gsub("\\([^d][^i][^s][^t][^=]*\\)", "", entry)
      scores <- sapply(symbols, function(symbol) {
        if (any(matches_criteria(symbol, cleaned_entry))) {
          score <- Open_Target_data$globalScore[Open_Target_data$symbol == symbol]
          if (length(score) > 0) matched_symbols <<- c(matched_symbols, symbol)
          return(score)
        }
        return(0)
      })
      scores <- unlist(scores)
      global_score_vector[i] <- if (length(scores) == 0) 0 else max(scores, na.rm = TRUE)
    } else {
      global_score_vector[i] <- 0
    }
  }
  cat("Number of variants with non-zero patho scores:", sum(global_score_vector > 0), "\n")
  cat("Number of unique gene symbols matched:", length(unique(matched_symbols)), "\n")

  pathogenic_indices <- which(global_score_vector > 0 & sapply(clnsig, function(x) !is.na(x) && x != "" && any(unlist(strsplit(x, "\\\\")) %in% c("Pathogenic", "Likely_pathogenic"))))

  # Tier definitions
  tier2_inc <- c("frameshift insertion", "frameshift deletion", "stopgain")
  gencode_exonic_info_tier2 <- get_exonic_indices(aGDS, "annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category", tier2_inc)
  refseq_exonic_info_tier2 <- get_exonic_indices(aGDS, "annotation/info/FunctionalAnnotation/refseq_exonic_category", tier2_inc)
  ucsc_exonic_info_tier2 <- get_exonic_indices(aGDS, "annotation/info/FunctionalAnnotation/ucsc_exonic_category", tier2_inc)

  tier3_inc <- c("nonsynonymous SNV", "nonframeshift deletion", "nonframeshift insertion", "stoploss")
  gencode_exonic_info_tier3 <- get_exonic_indices(aGDS, "annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category", tier3_inc)
  refseq_exonic_info_tier3 <- get_exonic_indices(aGDS, "annotation/info/FunctionalAnnotation/refseq_exonic_category", tier3_inc)
  ucsc_exonic_info_tier3 <- get_exonic_indices(aGDS, "annotation/info/FunctionalAnnotation/ucsc_exonic_category", tier3_inc)

  apc_protein_values <- seqGetData(aGDS, "annotation/info/FunctionalAnnotation/apc_protein_function_v3")
  apc_protein_values[is.na(apc_protein_values)] <- 0
  apc_protein_values <- pmin(apc_protein_values, 40)
  normalized_apc_protein_values <- (apc_protein_values - 0) / 40

  patho_score <- numeric(length(global_score_vector))
  patho_score_calc <- character(length(global_score_vector))
  valid_indices <- which(global_score_vector != 0)
  tier1_count <- 0; tier2_count <- 0; tier3_count <- 0; tier4_count <- 0

  for (i in seq_along(valid_indices)) {
    variant <- valid_indices[i]
    if (variant %in% pathogenic_indices) {
      patho_score[variant] = 80 + 20 * global_score_vector[variant]
      patho_score_calc[variant] = paste("Tier 1, 80 + 20 *", global_score_vector[variant])
      tier1_count <- tier1_count + 1
    } else if (variant %in% gencode_exonic_info_tier2 || variant %in% refseq_exonic_info_tier2 || variant %in% ucsc_exonic_info_tier2) {
      patho_score[variant] = 60 + 40 * global_score_vector[variant]
      patho_score_calc[variant] = paste("Tier 2, 60 + 40 *", global_score_vector[variant])
      tier2_count <- tier2_count + 1
    } else {
      score_tier3 <- 0; score_tier4 <- 0; calc_tier3 <- ""; calc_tier4 <- ""
      if (variant %in% gencode_exonic_info_tier3 || variant %in% refseq_exonic_info_tier3 || variant %in% ucsc_exonic_info_tier3) {
        score_tier3 = 20 + 80 * global_score_vector[variant]
        calc_tier3 = paste("Tier 3 = 20 + 80 *", global_score_vector[variant])
      }
      if (variant %in% valid_indices && apc_protein_values[which(valid_indices == variant)] > 1) {
        score_tier4 = 100 * ((0.5 * normalized_apc_protein_values[variant] + 0.5 * global_score_vector[variant]))
        calc_tier4 = paste("Tier 4 = 100 * ((0.5 *", normalized_apc_protein_values[variant], "+ 0.5 *", global_score_vector[variant], "))")
      }
      score_tier3 <- ifelse(is.na(score_tier3), 0, score_tier3)
      score_tier4 <- ifelse(is.na(score_tier4), 0, score_tier4)
      if (score_tier3 > score_tier4) {
        patho_score[variant] = score_tier3
        patho_score_calc[variant] = calc_tier3
        tier3_count <- tier3_count + 1
      } else {
        patho_score[variant] = score_tier4
        patho_score_calc[variant] = calc_tier4
        tier4_count <- tier4_count + 1
      }
    }
  }
  cat("Number of non-zero patho_score values:", sum(patho_score != 0), "\n")
  cat("Tier 1 variants:", tier1_count, "\n")
  cat("Tier 2 variants:", tier2_count, "\n")
  cat("Tier 3 variants:", tier3_count, "\n")
  cat("Tier 4 variants:", tier4_count, "\n")

  seqAddValue(aGDS, "annotation/info/patho_score", patho_score, replace = TRUE)
  seqAddValue(aGDS, "annotation/info/patho_score_calc", patho_score_calc, replace = TRUE)
  seqResetFilter(aGDS)

  sample_ids <- seqGetData(aGDS, "sample.id")
  for (sample_id in sample_ids) {
    # Filter to the current sample
    seqSetFilter(aGDS, sample.id = sample_id)

    # Get genotype matrix for the current sample
    genotypes <- seqGetData(aGDS, "genotype")

    # Determine valid variants for this sample
    valid_variants <- which(!is.na(genotypes[1, 1, ]) | !is.na(genotypes[2, 1, ]))

    # Get patho scores and filter to non-zero ones
    patho_scores <- seqGetData(aGDS, "annotation/info/patho_score")
    all_variants <- which(patho_scores > 0)
    all_variant_ids <- intersect(valid_variants, all_variants)

    # Apply final filter
    seqSetFilter(aGDS, sample.id = sample_id, variant.id = all_variant_ids)

    # Export to sample-specific chromosome file
    all_gdsfile <- paste0(sample_id, "_chr", chr, ".gds")
    seqExport(aGDS, all_gdsfile)

    # Reset filter
    seqResetFilter(aGDS)
  }

  seqClose(aGDS)
}

# Main script
main <- function() {
  cat("Main function started...\n")
  gds_files <- list.files(pattern = paste0("^", prefix, ".*\\.gds$"))
  chr_list <- sapply(gds_files, function(f) {
    matches <- regexec(pattern, f)
    match <- regmatches(f, matches)[[1]]
    if (length(match) > 1) return(match[2]) else return(NA)
  })
  names(gds_files) <- chr_list

  for (chr in chr_list[!is.na(chr_list)]) {
    gdsfile <- gds_files[[chr]]
    print(paste("Processing", gdsfile))
    new_gdsfile <- paste0(outprefix, chr, ".gds")
    file.copy(gdsfile, new_gdsfile, overwrite = TRUE)
    aGDS <- seqOpen(new_gdsfile, readonly = FALSE)
    seqResetFilter(aGDS)
    score_variants(aGDS, Open_Target_data, outprefix, chr, gdsfile)
  }

  # Ensure output directory exists
  dir.create("out", showWarnings = FALSE)

  # Move all final GDS files to out/
  file.copy(list.files(pattern = "_SNV_IMPACT\\.gds$"), "out", overwrite = TRUE)
  file.copy(list.files(pattern = "_chr[0-9XYM]+\\.gds$"), "out", overwrite = TRUE)
}

main()

merge_sample_gds_files <- function() {
  library(SeqArray)

  # Ensure output directory exists
  dir.create("out", showWarnings = FALSE)

  # List all GDS files that match the expected chromosome pattern
  gds_files <- list.files(pattern = ".*_chr[0-9XYM]+\\.gds$")

  # Extract sample IDs by removing the chromosome suffix
  sample_ids <- unique(sub("_chr[0-9XYM]+\\.gds$", "", gds_files))

  # Exclude generic or pipeline-generated prefixes
  sample_ids <- sample_ids[!sample_ids %in% c("merged", "anno_merged")]

  merge_sample_files <- function(sample_id) {
    sample_files <- list.files(pattern = paste0("^", sample_id, "_chr[0-9XYM]+\\.gds$"))
    if (length(sample_files) > 0) {
      merged_file <- paste0(sample_id, "_SNV_IMPACT.gds")
      message("Merging files for sample: ", sample_id)
      seqMerge(sample_files, merged_file, verbose = TRUE)

      # Move merged file to output directory
      file.copy(merged_file, file.path("out", merged_file), overwrite = TRUE)
    } else {
      message("No files found for sample: ", sample_id)
    }
  }

  for (sample_id in sample_ids) {
    merge_sample_files(sample_id)
  }
}

# Run the function
merge_sample_gds_files()
