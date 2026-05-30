#!/usr/bin/env Rscript

## Convert a flattened IMPACT-SNV/FAVOR parquet file into a SeqArray GDS.

if (requireNamespace("renv", quietly = TRUE)) {
  tryCatch(renv::activate(), error = function(e) message("renv activation failed: ", e$message))
}

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) stop("Package 'optparse' is required.")
  if (!requireNamespace("arrow", quietly = TRUE)) stop("Package 'arrow' is required.")
  if (!requireNamespace("SeqArray", quietly = TRUE)) stop("Bioconductor package 'SeqArray' is required.")
  if (!requireNamespace("gdsfmt", quietly = TRUE)) stop("Bioconductor package 'gdsfmt' is required.")
  library(optparse)
  library(arrow)
  library(SeqArray)
  library(gdsfmt)
})

option_list <- list(
  make_option(c("-i", "--input"), dest = "input", type = "character", help = "Input flat parquet"),
  make_option(c("-o", "--output"), dest = "output", type = "character", help = "Output SeqArray GDS path"),
  make_option(c("-s", "--sample-id"), dest = "sample_id", type = "character", default = NULL, help = "Sample ID"),
  make_option(c("--force"), dest = "force", action = "store_true", default = FALSE, help = "Overwrite output GDS"),
  make_option(c("--keep-temp-vcf"), dest = "keep_temp_vcf", action = "store_true", default = FALSE, help = "Keep temporary VCF"),
  make_option(c("--temp-vcf"), dest = "temp_vcf", type = "character", default = NULL, help = "Optional temporary VCF path"),
  make_option(c("--no-optimize"), dest = "no_optimize", action = "store_true", default = FALSE, help = "Skip seqOptimize"),
  make_option(c("--storage-option"), dest = "storage_option", type = "character", default = "ZIP_RA", help = "SeqArray storage option"),
  make_option(c("--verbose"), dest = "verbose", action = "store_true", default = FALSE, help = "Verbose output")
)
opt <- parse_args(OptionParser(option_list = option_list))

required_columns <- c(
  "sample_id", "variant_id", "chromosome", "position", "ref", "alt", "allele",
  "dosage", "maf", "vid", "variant_vcf", "matched_gene", "matched_gene_score",
  "genecode_comprehensive_info", "genecode_comprehensive_exonic_category",
  "refseq_exonic_category", "ucsc_exonic_category", "clnsig",
  "apc_protein_function_v3"
)

fail <- function(...) stop(paste0(...), call. = FALSE)
log_msg <- function(...) cat(paste0(..., "\n"))
as_chr <- function(x, default = "") { y <- as.character(x); y[is.na(y)] <- default; y }
as_num <- function(x) suppressWarnings(as.numeric(x))
as_int <- function(x) suppressWarnings(as.integer(x))
normalize_chrom <- function(x) { y <- as_chr(x); y <- sub("^chr", "", y, ignore.case = TRUE); y <- sub("\\.0$", "", y); y }
vcf_scalar <- function(x, default = ".") { y <- as_chr(x, default); y[y == ""] <- default; y }
vcf_filter <- function(x) { y <- as_chr(x, "PASS"); y[y == ""] <- "PASS"; y }
vcf_qual <- function(x) { y <- as_chr(x, "."); y[y == ""] <- "."; y }

dosage_to_gt_string <- function(dosage) {
  d <- as.integer(round(dosage))
  if (is.na(d)) return("./.")
  if (d <= 0L) return("0/0")
  if (d == 1L) return("0/1")
  "1/1"
}

ensure_parent_dir <- function(path) {
  parent <- dirname(normalizePath(path, mustWork = FALSE))
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
}

validate_input <- function(df) {
  missing <- setdiff(required_columns, names(df))
  if (length(missing) > 0) fail("Input flat file is missing required columns: ", paste(missing, collapse = ", "))
  if (nrow(df) == 0) fail("Input flat file has zero rows. Nothing to write to GDS.")
  if (any(is.na(df$position))) fail("Input contains NA positions.")
  if (any(is.na(df$dosage))) fail("Input contains NA dosage values.")
}

resolve_storage_option <- function(storage_option_name) {
  if (is.null(storage_option_name) || length(storage_option_name) == 0L || is.na(storage_option_name) || storage_option_name == "") {
    storage_option_name <- "ZIP_RA"
  }
  obj <- tryCatch(
    SeqArray::seqStorageOption(as.character(storage_option_name)),
    error = function(e) fail("Could not construct SeqArray storage option '", storage_option_name, "': ", e$message)
  )
  if (!inherits(obj, "SeqGDSStorageClass")) {
    fail("SeqArray::seqStorageOption('", storage_option_name, "') did not return a SeqGDSStorageClass object.")
  }
  obj
}

add_info <- function(gds, name, value, replace = TRUE) {
  SeqArray::seqAddValue(gds, paste0("annotation/info/", name), value, replace = replace, verbose = FALSE)
}

add_func <- function(gds, name, value, replace = TRUE) {
  SeqArray::seqAddValue(gds, paste0("annotation/info/FunctionalAnnotation/", name), value, replace = replace, verbose = FALSE)
}

ensure_functional_annotation_folder <- function(gds) {
  info_node <- gdsfmt::index.gdsn(gds, "annotation/info", silent = TRUE)
  if (is.null(info_node)) {
    ann_node <- gdsfmt::index.gdsn(gds, "annotation", silent = TRUE)
    if (is.null(ann_node)) ann_node <- gdsfmt::addfolder.gdsn(gds, "annotation")
    info_node <- gdsfmt::addfolder.gdsn(ann_node, "info")
  }
  func_node <- gdsfmt::index.gdsn(info_node, "FunctionalAnnotation", silent = TRUE)
  if (is.null(func_node)) gdsfmt::addfolder.gdsn(info_node, "FunctionalAnnotation")
  invisible(TRUE)
}

write_temp_vcf <- function(df, sample_id, temp_vcf) {
  log_msg("Writing temporary VCF: ", temp_vcf)
  ensure_parent_dir(temp_vcf)
  nvar <- nrow(df)
  id <- if ("rsid" %in% names(df)) vcf_scalar(df$rsid, ".") else rep(".", nvar)
  qual <- if ("qual" %in% names(df)) vcf_qual(df$qual) else rep(".", nvar)
  filter <- if ("filter" %in% names(df)) vcf_filter(df$filter) else rep("PASS", nvar)
  gt <- vapply(df$dosage, dosage_to_gt_string, character(1L))
  ds <- as_chr(df$dosage)
  sample_value <- paste0(gt, ":", ds)

  con <- file(temp_vcf, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("##fileformat=VCFv4.2", con)
  writeLines("##source=IMPACT_favor_flat_to_seqarray_gds", con)
  writeLines('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype reconstructed from IMPACT-SNV dosage">', con)
  writeLines('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage from IMPACT-SNV genotype parquet">', con)
  for (chrom in unique(as_chr(df$chromosome))) writeLines(paste0("##contig=<ID=", chrom, ">"), con)
  writeLines(paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_id), collapse = "\t"), con)
  lines <- paste(as_chr(df$chromosome), as_int(df$position), id, as_chr(df$ref), as_chr(df$alt), qual, filter, ".", "GT:DS", sample_value, sep = "\t")
  writeLines(lines, con)
  invisible(temp_vcf)
}

append_annotations <- function(gds_path, df, input_path) {
  log_msg("Appending IMPACT/FAVOR annotations to GDS...")
  gds <- SeqArray::seqOpen(gds_path, readonly = FALSE)
  on.exit(try(SeqArray::seqClose(gds), silent = TRUE), add = TRUE)
  ensure_functional_annotation_folder(gds)
  nvar <- nrow(df)

  # Use sequential integer variant IDs. This is safer than trusting legacy/string IDs.
  SeqArray::seqAddValue(gds, "variant.id", seq_len(nvar), replace = TRUE, verbose = FALSE)
  add_info(gds, "variant_vcf", as_chr(df$variant_vcf))
  add_info(gds, "vid", as_chr(df$vid))
  add_info(gds, "ref_vcf", as_chr(df$ref))
  add_info(gds, "alt_vcf", as_chr(df$alt))
  add_info(gds, "maf", as_num(df$maf))
  add_info(gds, "dosage", as_num(df$dosage))
  add_info(gds, "matched_gene", as_chr(df$matched_gene))
  add_info(gds, "matched_gene_score", as_num(df$matched_gene_score))
  if ("matched_gene_all" %in% names(df)) add_info(gds, "matched_gene_all", as_chr(df$matched_gene_all))
  if ("matched_gene_score_all" %in% names(df)) add_info(gds, "matched_gene_score_all", as_chr(df$matched_gene_score_all))
  if ("matched_gene_source" %in% names(df)) add_info(gds, "matched_gene_source", as_chr(df$matched_gene_source))

  add_func(gds, "clnsig", as_chr(df$clnsig))
  add_func(gds, "genecode_comprehensive_info", as_chr(df$genecode_comprehensive_info))
  add_func(gds, "genecode_comprehensive_exonic_category", as_chr(df$genecode_comprehensive_exonic_category))
  add_func(gds, "refseq_exonic_category", as_chr(df$refseq_exonic_category))
  add_func(gds, "ucsc_exonic_category", as_chr(df$ucsc_exonic_category))
  add_func(gds, "apc_protein_function_v3", as_num(df$apc_protein_function_v3))

  extra_map <- list(
    matched_gene = "matched_gene",
    matched_gene_score = "matched_gene_score",
    gencode_genes = "gencode_genes",
    gencode_region_type = "gencode_region_type",
    gencode_consequence = "gencode_consequence",
    refseq_consequence = "refseq_consequence",
    ucsc_consequence = "ucsc_consequence",
    clndn = "clndn",
    clinvar_gene = "clinvar_gene",
    spliceai_max_ds = "spliceai_max_ds",
    alphamissense_max_pathogenicity = "alphamissense_max_pathogenicity",
    revel = "revel",
    gnomad_genome_af = "gnomad_genome_af",
    gnomad_exome_af = "gnomad_exome_af",
    bravo_af = "bravo_af",
    tg_all = "tg_all",
    cadd_phred = "cadd_phred"
  )
  numeric_extras <- c("matched_gene_score", "spliceai_max_ds", "alphamissense_max_pathogenicity", "revel", "gnomad_genome_af", "gnomad_exome_af", "bravo_af", "tg_all", "cadd_phred")
  for (node_name in names(extra_map)) {
    col <- extra_map[[node_name]]
    if (col %in% names(df)) {
      value <- if (col %in% numeric_extras) as_num(df[[col]]) else as_chr(df[[col]])
      add_func(gds, node_name, value)
    }
  }
  add_info(gds, "source_flat_file", rep(normalizePath(input_path, mustWork = FALSE), nvar))
  add_info(gds, "gds_contract_version", rep("impact_snv_v1_flat_contract", nvar))
  invisible(TRUE)
}

create_gds <- function(input_path, output_path, sample_id = NULL, force = FALSE, keep_temp_vcf = FALSE, temp_vcf = NULL, optimize = TRUE, storage_option = "ZIP_RA", verbose = FALSE) {
  if (!file.exists(input_path)) fail("Input parquet not found: ", input_path)
  if (file.exists(output_path)) {
    if (!force) fail("Output exists; use --force to overwrite: ", output_path)
    unlink(output_path, force = TRUE)
  }
  ensure_parent_dir(output_path)
  log_msg("Resolved options:")
  log_msg("  input:           ", input_path)
  log_msg("  output:          ", output_path)
  log_msg("  sample_id:       ", ifelse(is.null(sample_id), "", sample_id))
  log_msg("  keep_temp_vcf:   ", keep_temp_vcf)
  log_msg("  temp_vcf:        ", ifelse(is.null(temp_vcf) || temp_vcf == "", "", temp_vcf))
  log_msg("  storage_option:  ", storage_option)
  log_msg("  optimize:        ", optimize)

  log_msg("Reading flat parquet: ", input_path)
  df <- as.data.frame(arrow::read_parquet(input_path))
  validate_input(df)
  if (is.null(sample_id) || sample_id == "") {
    sample_values <- unique(as_chr(df$sample_id))
    if (length(sample_values) != 1L) fail("--sample-id was not provided and input contains multiple sample_id values: ", paste(sample_values, collapse = ", "))
    sample_id <- sample_values[[1L]]
  }
  df <- df[as_chr(df$sample_id) == sample_id, , drop = FALSE]
  if (nrow(df) == 0) fail("No rows remain for sample_id: ", sample_id)
  df$chromosome <- normalize_chrom(df$chromosome)
  df$position <- as_int(df$position)
  df <- df[order(df$chromosome, df$position, as_chr(df$ref), as_chr(df$alt)), , drop = FALSE]
  nvar <- nrow(df)
  log_msg("Preparing SeqArray GDS")
  log_msg("  sample_id: ", sample_id)
  log_msg("  variants:  ", nvar)

  if (is.null(temp_vcf) || temp_vcf == "") temp_vcf <- tempfile(pattern = paste0(sample_id, "_"), fileext = ".vcf")
  write_temp_vcf(df, sample_id, temp_vcf)
  storage_obj <- resolve_storage_option(storage_option)
  log_msg("Using SeqArray storage option: ", storage_option)
  log_msg("Converting temporary VCF to SeqArray GDS: ", output_path)
  SeqArray::seqVCF2GDS(vcf.fn = temp_vcf, out.fn = output_path, storage.option = storage_obj, fmt.import = c("GT", "DS"), optimize = FALSE, verbose = TRUE)
  append_annotations(output_path, df, input_path)
  if (optimize) {
    log_msg("Optimizing GDS...")
    SeqArray::seqOptimize(output_path)
  }
  if (!keep_temp_vcf) unlink(temp_vcf, force = TRUE) else log_msg("Temporary VCF retained: ", temp_vcf)

  log_msg("Running post-write sanity check...")
  g <- SeqArray::seqOpen(output_path)
  on.exit(try(SeqArray::seqClose(g), silent = TRUE), add = TRUE)
  observed_sample <- SeqArray::seqGetData(g, "sample.id")
  observed_variants <- length(SeqArray::seqGetData(g, "variant.id"))
  observed_dosage <- SeqArray::seqGetData(g, "annotation/info/dosage")
  observed_gene <- SeqArray::seqGetData(g, "annotation/info/matched_gene")
  if (!identical(as.character(observed_sample), as.character(sample_id))) fail("Post-write sample.id check failed")
  if (observed_variants != nvar) fail("Post-write variant count check failed: expected ", nvar, ", observed ", observed_variants)
  if (length(observed_dosage) != nvar || length(observed_gene) != nvar) fail("Post-write annotation length check failed")
  SeqArray::seqClose(g)
  log_msg("Done: ", output_path)
  invisible(output_path)
}

if (is.null(opt$input) || is.null(opt$output)) fail("Both --input and --output are required. Use --help for usage.")
create_gds(input_path = opt$input, output_path = opt$output, sample_id = opt$sample_id, force = opt$force, keep_temp_vcf = opt$keep_temp_vcf, temp_vcf = opt$temp_vcf, optimize = !opt$no_optimize, storage_option = opt$storage_option, verbose = opt$verbose)
