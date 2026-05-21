#!/usr/bin/env Rscript
# add_impact_vis_compat_nodes.R
#
# Add IMPACT-VIS compatibility aliases/placeholders to a finalized IMPACT-SNV GDS.
#
# Critical IMPACT-VIS compatibility node:
#   annotation/info/FunctionalAnnotation/VarInfo
#
# IMPORTANT:
#   VarInfo is the canonical variant key in chr-pos-ref-alt format, e.g.
#   8-9719889-T-A
#
# This script derives VarInfo from chromosome, position, REF, and ALT. REF/ALT are
# taken from annotation/info/ref_vcf and annotation/info/alt_vcf if available;
# otherwise they are parsed from the SeqArray allele node, which is usually
# stored as REF,ALT.

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) stop("Package 'optparse' is required")
  if (!requireNamespace("SeqArray", quietly = TRUE)) stop("Bioconductor package 'SeqArray' is required")
  if (!requireNamespace("gdsfmt", quietly = TRUE)) stop("Bioconductor package 'gdsfmt' is required")
  library(optparse)
  library(SeqArray)
  library(gdsfmt)
})

option_list <- list(
  make_option(c("-i", "--input-gds"), dest = "input_gds", type = "character",
              help = "Input finalized *_SNV_IMPACT.gds"),
  make_option(c("-o", "--output-gds"), dest = "output_gds", type = "character", default = NULL,
              help = "Optional output GDS. If omitted with --in-place, modifies input in place."),
  make_option(c("--in-place"), dest = "in_place", action = "store_true", default = FALSE,
              help = "Modify input GDS in place"),
  make_option(c("--force"), dest = "force", action = "store_true", default = FALSE,
              help = "Overwrite output GDS if it exists"),
  make_option(c("--no-optimize"), dest = "no_optimize", action = "store_true", default = FALSE,
              help = "Skip seqOptimize() after adding nodes")
)
opt <- parse_args(OptionParser(option_list = option_list))

fail <- function(...) stop(paste0(...), call. = FALSE)
msg <- function(...) cat(paste0(..., "\n"))

as_chr <- function(x, default = "") {
  y <- as.character(x)
  y[is.na(y)] <- default
  y
}

get_or_default <- function(gds, node, default) {
  tryCatch(SeqArray::seqGetData(gds, node), error = function(e) default)
}

node_exists <- function(gds, node) {
  !is.null(gdsfmt::index.gdsn(gds, node, silent = TRUE))
}

ensure_folder <- function(gds, path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  cur <- gds
  current_path <- character(0)
  for (p in parts) {
    current_path <- c(current_path, p)
    existing <- gdsfmt::index.gdsn(gds, paste(current_path, collapse = "/"), silent = TRUE)
    if (is.null(existing)) {
      cur <- gdsfmt::addfolder.gdsn(cur, p)
    } else {
      cur <- existing
    }
  }
  invisible(cur)
}

add_seq_value <- function(gds, node, value) {
  SeqArray::seqAddValue(gds, node, value, replace = TRUE, verbose = FALSE)
}

parse_ref_alt_from_allele <- function(allele) {
  allele <- as_chr(allele)
  split <- strsplit(allele, ",", fixed = TRUE)
  ref <- vapply(split, function(x) if (length(x) >= 1) x[[1]] else "", character(1))
  alt <- vapply(split, function(x) if (length(x) >= 2) paste(x[-1], collapse = ",") else "", character(1))
  list(ref = ref, alt = alt)
}

make_varinfo <- function(gds, n) {
  chrom <- as_chr(SeqArray::seqGetData(gds, "chromosome"))
  chrom <- sub("^chr", "", chrom, ignore.case = TRUE)
  pos <- as_chr(SeqArray::seqGetData(gds, "position"))

  # Prefer explicit VCF REF/ALT annotations written by favor_flat_to_seqarray_gds.R.
  ref <- as_chr(get_or_default(gds, "annotation/info/ref_vcf", rep("", n)))
  alt <- as_chr(get_or_default(gds, "annotation/info/alt_vcf", rep("", n)))

  # Fallback to SeqArray allele node, usually REF,ALT.
  missing_ref_alt <- ref == "" | alt == ""
  if (any(missing_ref_alt)) {
    parsed <- parse_ref_alt_from_allele(SeqArray::seqGetData(gds, "allele"))
    ref[missing_ref_alt & ref == ""] <- parsed$ref[missing_ref_alt & ref == ""]
    alt[missing_ref_alt & alt == ""] <- parsed$alt[missing_ref_alt & alt == ""]
  }

  if (any(chrom == "" | pos == "" | ref == "" | alt == "")) {
    bad <- which(chrom == "" | pos == "" | ref == "" | alt == "")
    fail("Could not construct VarInfo for ", length(bad), " variants; first bad index: ", bad[[1]])
  }

  paste(chrom, pos, ref, alt, sep = "-")
}

main <- function() {
  if (is.null(opt$input_gds)) fail("--input-gds is required")
  if (!file.exists(opt$input_gds)) fail("Input GDS not found: ", opt$input_gds)

  if (opt$in_place) {
    target <- opt$input_gds
  } else {
    if (is.null(opt$output_gds)) fail("Provide --output-gds or use --in-place")
    target <- opt$output_gds
    if (file.exists(target)) {
      if (!opt$force) fail("Output exists; use --force to overwrite: ", target)
      unlink(target, force = TRUE)
    }
    dir.create(dirname(normalizePath(target, mustWork = FALSE)), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(opt$input_gds, target, overwrite = TRUE)) fail("Could not copy input to output: ", target)
  }

  msg("Adding IMPACT-VIS compatibility nodes to: ", target)
  gds <- SeqArray::seqOpen(target, readonly = FALSE)
  on.exit(try(SeqArray::seqClose(gds), silent = TRUE), add = TRUE)

  n <- length(SeqArray::seqGetData(gds, "variant.id"))
  ensure_folder(gds, "annotation/info/FunctionalAnnotation")

  # Critical IMPACT-VIS compatibility node: chr-pos-ref-alt.
  varinfo <- make_varinfo(gds, n)
  add_seq_value(gds, "annotation/info/FunctionalAnnotation/VarInfo", varinfo)

  # Ensure several commonly expected legacy FunctionalAnnotation nodes exist.
  # These are safe placeholders if the FAVOR-flat source did not provide them.
  if (!node_exists(gds, "annotation/info/FunctionalAnnotation/clnsig")) {
    add_seq_value(gds, "annotation/info/FunctionalAnnotation/clnsig", rep("", n))
  }
  if (!node_exists(gds, "annotation/info/FunctionalAnnotation/bravo_af")) {
    add_seq_value(gds, "annotation/info/FunctionalAnnotation/bravo_af", rep(NA_real_, n))
  }
  if (!node_exists(gds, "annotation/info/FunctionalAnnotation/aloft_prediction")) {
    add_seq_value(gds, "annotation/info/FunctionalAnnotation/aloft_prediction", rep("", n))
  }

  required <- c(
    "annotation/info/impact_score",
    "annotation/info/impact_score_calc",
    "annotation/info/tier",
    "annotation/info/FunctionalAnnotation/VarInfo"
  )
  for (node in required) {
    x <- SeqArray::seqGetData(gds, node)
    if (length(x) != n) fail("Node ", node, " has length ", length(x), ", expected ", n)
  }

  msg("VarInfo preview:")
  print(utils::head(varinfo, 10))

  SeqArray::seqClose(gds)
  gds <- NULL

  if (!opt$no_optimize) {
    msg("Optimizing GDS...")
    SeqArray::seqOptimize(target)
  }

  msg("Done: ", target)
}

main()
