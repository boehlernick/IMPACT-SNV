#!/usr/bin/env Rscript
# impact_prioritize_gds.R
#
# Add IMPACT scores, tier labels, and ClinVar significance flags to a
# pre-prioritization IMPACT/FAVOR SeqArray GDS file.
#
# Input:
#   A per-sample GDS produced by favor_flat_to_seqarray_gds.R / build_favor_sample_gds.py
#   and a GeneList.txt file with columns: symbol, globalScore.
#
# Output:
#   A final per-sample *_SNV_IMPACT.gds file ready for IMPACT-VIS-style
#   downstream use.
#
# Deprecated patho_score/patho_score_calc nodes are intentionally NOT generated.
# All scoring logic writes directly to impact_score/impact_score_calc.

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) stop("Package 'optparse' is required.")
  if (!requireNamespace("SeqArray", quietly = TRUE)) stop("Bioconductor package 'SeqArray' is required.")
  if (!requireNamespace("gdsfmt", quietly = TRUE)) stop("Bioconductor package 'gdsfmt' is required.")
  if (!requireNamespace("stringi", quietly = TRUE)) stop("Package 'stringi' is required.")
  library(optparse)
  library(SeqArray)
  library(gdsfmt)
  library(stringi)
})

option_list <- list(
  make_option(c("-i", "--input-gds"), dest = "input_gds", type = "character", help = "Input pre-prioritization GDS"),
  make_option(c("-g", "--gene-list"), dest = "gene_list", type = "character", help = "GeneList.txt with symbol and globalScore columns"),
  make_option(c("-o", "--output-gds"), dest = "output_gds", type = "character", help = "Output finalized *_SNV_IMPACT.gds"),
  make_option(c("--force"), dest = "force", action = "store_true", default = FALSE, help = "Overwrite output GDS if it exists"),
  make_option(c("--no-optimize"), dest = "no_optimize", action = "store_true", default = FALSE, help = "Skip seqOptimize() after adding scoring nodes"),
  make_option(c("--verbose"), dest = "verbose", action = "store_true", default = FALSE, help = "Print additional information")
)
opt <- parse_args(OptionParser(option_list = option_list))

fail <- function(...) stop(paste0(...), call. = FALSE)
msg <- function(...) cat(paste0(..., "\n"))

as_chr <- function(x, default = "") {
  y <- as.character(x)
  y[is.na(y)] <- default
  y
}

read_genelist <- function(path) {
  if (!file.exists(path)) fail("Gene list not found: ", path)
  df <- read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  required <- c("symbol", "globalScore")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) fail("Gene list missing required columns: ", paste(missing, collapse = ", "))
  df$symbol <- as.character(df$symbol)
  df$globalScore <- suppressWarnings(as.numeric(df$globalScore))
  df <- df[!is.na(df$symbol) & df$symbol != "", , drop = FALSE]
  if (nrow(df) == 0) fail("Gene list contains no usable symbols: ", path)
  score <- tapply(df$globalScore, df$symbol, function(x) max(x, na.rm = TRUE))
  score[is.infinite(score)] <- NA_real_
  score
}

split_gene_string <- function(x) {
  x <- as_chr(x)
  pieces <- unlist(strsplit(x, "[|,;]", perl = TRUE), use.names = FALSE)
  pieces <- trimws(pieces)
  pieces[pieces != ""]
}

score_from_genes <- function(matched_gene_all, matched_gene, gene_scores) {
  n <- length(matched_gene)
  out_score <- rep(0, n)
  out_gene <- rep("", n)
  for (i in seq_len(n)) {
    genes <- character(0)
    if (!is.na(matched_gene_all[i]) && matched_gene_all[i] != "") genes <- split_gene_string(matched_gene_all[i])
    if (length(genes) == 0 && !is.na(matched_gene[i]) && matched_gene[i] != "") genes <- matched_gene[i]
    genes <- genes[genes %in% names(gene_scores)]
    if (length(genes) > 0) {
      scores <- as.numeric(gene_scores[genes])
      scores[is.na(scores)] <- 0
      j <- which.max(scores)
      out_score[i] <- scores[j]
      out_gene[i] <- genes[j]
    }
  }
  list(score = out_score, gene = out_gene)
}

clinvar_tokens <- function(x) {
  x <- as_chr(x)
  x <- gsub("_", " ", x)
  toks <- unlist(strsplit(x, "\\\\|\\||;|,", perl = TRUE), use.names = FALSE)
  toks <- tolower(trimws(toks))
  toks[toks != "" & toks != "na"]
}

has_clinvar <- function(clnsig, labels) {
  labels <- tolower(labels)
  vapply(clnsig, function(x) {
    toks <- clinvar_tokens(x)
    any(toks %in% labels)
  }, logical(1L))
}

category_any <- function(categories, include) {
  include_l <- tolower(include)
  Reduce(`|`, lapply(categories, function(x) tolower(as_chr(x)) %in% include_l))
}

ensure_info_folder <- function(gds) {
  ann <- gdsfmt::index.gdsn(gds, "annotation", silent = TRUE)
  if (is.null(ann)) ann <- gdsfmt::addfolder.gdsn(gds, "annotation")
  info <- gdsfmt::index.gdsn(ann, "info", silent = TRUE)
  if (is.null(info)) info <- gdsfmt::addfolder.gdsn(ann, "info")
  info
}

add_clnsig_flags <- function(gds, clnsig) {
  info <- ensure_info_folder(gds)
  flags_node <- gdsfmt::index.gdsn(info, "clnsig_flags", silent = TRUE)
  if (is.null(flags_node)) flags_node <- gdsfmt::addfolder.gdsn(info, "clnsig_flags")

  label_to_flag <- c(
    "Pathogenic" = "pathogenic",
    "Likely pathogenic" = "likely_pathogenic",
    "Uncertain significance" = "uncertain_significance",
    "Likely benign" = "likely_benign",
    "Benign" = "benign",
    "Pathogenic, low penetrance" = "pathogenic_low_penetrance",
    "Likely pathogenic, low penetrance" = "likely_pathogenic_low_penetrance",
    "Established risk allele" = "established_risk_allele",
    "Likely risk allele" = "likely_risk_allele",
    "Uncertain risk allele" = "uncertain_risk_allele",
    "affects" = "affects",
    "association" = "association",
    "drug response" = "drug_response",
    "confers sensitivity" = "confers_sensitivity",
    "protective" = "protective",
    "other" = "other",
    "Conflicting interpretations of pathogenicity" = "conflicting_interpretations_of_pathogenicity",
    "not provided" = "not_provided"
  )

  standard <- tolower(names(label_to_flag))
  token_list <- lapply(clnsig, clinvar_tokens)

  for (k in seq_along(label_to_flag)) {
    label <- tolower(names(label_to_flag)[k])
    node_name <- unname(label_to_flag[k])
    values <- vapply(token_list, function(toks) label %in% toks, logical(1L))
    gdsfmt::add.gdsn(flags_node, node_name, val = values, storage = "bit1", replace = TRUE)
  }

  other_values <- vapply(token_list, function(toks) {
    length(toks) > 0 && ("other" %in% toks || any(!(toks %in% standard)))
  }, logical(1L))
  gdsfmt::add.gdsn(flags_node, "other", val = other_values, storage = "bit1", replace = TRUE)
  invisible(TRUE)
}

score_variants <- function(gds, gene_scores, verbose = FALSE) {
  variant_id <- SeqArray::seqGetData(gds, "variant.id")
  n <- length(variant_id)

  get_node <- function(node, default = rep("", n)) {
    tryCatch(SeqArray::seqGetData(gds, node), error = function(e) default)
  }

  clnsig <- as_chr(get_node("annotation/info/FunctionalAnnotation/clnsig"))
  matched_gene <- as_chr(get_node("annotation/info/matched_gene"))
  matched_gene_all <- as_chr(get_node("annotation/info/matched_gene_all"))
  stored_gene_score <- suppressWarnings(as.numeric(get_node("annotation/info/matched_gene_score", rep(NA_real_, n))))

  score_info <- score_from_genes(matched_gene_all, matched_gene, gene_scores)
  gene_score <- score_info$score
  best_gene <- score_info$gene
  use_stored <- gene_score == 0 & !is.na(stored_gene_score)
  gene_score[use_stored] <- stored_gene_score[use_stored]
  best_gene[best_gene == "" & matched_gene != ""] <- matched_gene[best_gene == "" & matched_gene != ""]
  gene_score[is.na(gene_score)] <- 0

  gencode_cat <- as_chr(get_node("annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category"))
  refseq_cat <- as_chr(get_node("annotation/info/FunctionalAnnotation/refseq_exonic_category"))
  ucsc_cat <- as_chr(get_node("annotation/info/FunctionalAnnotation/ucsc_exonic_category"))
  apc <- suppressWarnings(as.numeric(get_node("annotation/info/FunctionalAnnotation/apc_protein_function_v3", rep(NA_real_, n))))
  apc[is.na(apc)] <- 0

  category_set <- list(gencode_cat, refseq_cat, ucsc_cat)
  tier1 <- has_clinvar(clnsig, c("pathogenic", "likely pathogenic")) & gene_score > 0
  tier2 <- category_any(category_set, c("frameshift insertion", "frameshift deletion", "stopgain")) & gene_score > 0
  tier3 <- category_any(category_set, c("nonsynonymous SNV", "nonframeshift deletion", "nonframeshift insertion", "stoploss")) & gene_score > 0
  tier4 <- apc > 1 & gene_score > 0
  apc_norm <- pmin(apc, 40) / 40

  impact_score <- rep(0, n)
  impact_calc <- rep("Tier 0, no Tier 1-4 criteria met", n)
  tier <- rep(0L, n)

  idx <- which(tier1)
  impact_score[idx] <- 80 + 20 * gene_score[idx]
  impact_calc[idx] <- paste0("Tier 1, 80 + 20 * ", signif(gene_score[idx], 6), " [", best_gene[idx], "]")
  tier[idx] <- 1L

  idx <- which(tier == 0L & tier2)
  impact_score[idx] <- 60 + 40 * gene_score[idx]
  impact_calc[idx] <- paste0("Tier 2, 60 + 40 * ", signif(gene_score[idx], 6), " [", best_gene[idx], "]")
  tier[idx] <- 2L

  remaining <- which(tier == 0L & (tier3 | tier4))
  if (length(remaining) > 0) {
    score3 <- rep(0, length(remaining))
    score4 <- rep(0, length(remaining))
    idx_t3 <- which(tier3[remaining])
    idx_t4 <- which(tier4[remaining])
    score3[idx_t3] <- 20 + 80 * gene_score[remaining[idx_t3]]
    score4[idx_t4] <- 100 * (0.5 * apc_norm[remaining[idx_t4]] + 0.5 * gene_score[remaining[idx_t4]])
    choose3 <- score3 >= score4
    idx3 <- remaining[choose3 & score3 > 0]
    idx4 <- remaining[!choose3 & score4 > 0]
    impact_score[idx3] <- 20 + 80 * gene_score[idx3]
    impact_calc[idx3] <- paste0("Tier 3, 20 + 80 * ", signif(gene_score[idx3], 6), " [", best_gene[idx3], "]")
    tier[idx3] <- 3L
    impact_score[idx4] <- 100 * (0.5 * apc_norm[idx4] + 0.5 * gene_score[idx4])
    impact_calc[idx4] <- paste0("Tier 4, 100 * ((0.5 * ", signif(apc_norm[idx4], 6), ") + (0.5 * ", signif(gene_score[idx4], 6), ")) [", best_gene[idx4], "]")
    tier[idx4] <- 4L
  }

  # Only non-deprecated score nodes are generated.
  SeqArray::seqAddValue(gds, "annotation/info/impact_score", impact_score, replace = TRUE, verbose = FALSE)
  SeqArray::seqAddValue(gds, "annotation/info/impact_score_calc", impact_calc, replace = TRUE, verbose = FALSE)
  SeqArray::seqAddValue(gds, "annotation/info/tier", as.integer(tier), replace = TRUE, verbose = FALSE)
  SeqArray::seqAddValue(gds, "annotation/info/scoring_gene", best_gene, replace = TRUE, verbose = FALSE)
  SeqArray::seqAddValue(gds, "annotation/info/scoring_gene_score", gene_score, replace = TRUE, verbose = FALSE)
  add_clnsig_flags(gds, clnsig)

  list(
    variant_count = n,
    nonzero_impact_score_count = sum(impact_score > 0, na.rm = TRUE),
    tier0 = sum(tier == 0L),
    tier1 = sum(tier == 1L),
    tier2 = sum(tier == 2L),
    tier3 = sum(tier == 3L),
    tier4 = sum(tier == 4L),
    clnsig_nonempty = sum(clnsig != ""),
    gene_score_positive = sum(gene_score > 0),
    max_impact_score = max(impact_score, na.rm = TRUE)
  )
}

copy_input <- function(input_gds, output_gds, force = FALSE) {
  if (!file.exists(input_gds)) fail("Input GDS not found: ", input_gds)
  if (file.exists(output_gds)) {
    if (!force) fail("Output GDS exists; use --force to overwrite: ", output_gds)
    unlink(output_gds, force = TRUE)
  }
  dir.create(dirname(normalizePath(output_gds, mustWork = FALSE)), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(input_gds, output_gds, overwrite = TRUE)
  if (!ok) fail("Could not copy input GDS to output: ", output_gds)
}

main <- function() {
  if (is.null(opt$input_gds) || is.null(opt$output_gds) || is.null(opt$gene_list)) {
    fail("--input-gds, --output-gds, and --gene-list are required")
  }
  msg("Input GDS:  ", opt$input_gds)
  msg("Output GDS: ", opt$output_gds)
  msg("Gene list:  ", opt$gene_list)
  gene_scores <- read_genelist(opt$gene_list)
  msg("Loaded ", length(gene_scores), " GeneList symbols")
  copy_input(opt$input_gds, opt$output_gds, force = opt$force)

  gds <- SeqArray::seqOpen(opt$output_gds, readonly = FALSE)
  on.exit(try(SeqArray::seqClose(gds), silent = TRUE), add = TRUE)
  summary <- score_variants(gds, gene_scores, verbose = opt$verbose)
  SeqArray::seqClose(gds)
  gds <- NULL

  if (!opt$no_optimize) {
    msg("Optimizing scored GDS...")
    SeqArray::seqOptimize(opt$output_gds)
  }
  msg("IMPACT prioritization complete")
  msg("  variants:                     ", summary$variant_count)
  msg("  nonzero impact_score variants: ", summary$nonzero_impact_score_count)
  msg("  tier0:                        ", summary$tier0)
  msg("  tier1:                        ", summary$tier1)
  msg("  tier2:                        ", summary$tier2)
  msg("  tier3:                        ", summary$tier3)
  msg("  tier4:                        ", summary$tier4)
  msg("  ClinVar non-empty:            ", summary$clnsig_nonempty)
  msg("  positive gene scores:         ", summary$gene_score_positive)
  msg("  max impact_score:             ", signif(summary$max_impact_score, 6))
  msg("Done: ", opt$output_gds)
}
main()
