#!/usr/bin/env Rscript

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i + 1 <= length(args)) args[[i + 1]] else ""
    if (startsWith(key, "--")) {
      parsed[[sub("^--", "", key)]] <- value
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  parsed
}

args <- parse_args()
input_path <- args[["input"]]
output_path <- args[["output"]]
reference_genome_build <- args[["reference-genome-build"]]
dry_run <- args[["dry-run"]]

if (is.null(input_path) || is.null(output_path)) {
  stop("Both --input and --output are required.")
}

table <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("variant.id", "chromosome", "position", "ref", "alt", "canonical_variant_key")
missing_columns <- setdiff(required_columns, names(table))

if (length(missing_columns) > 0) {
  stop(sprintf("Missing required column(s): %s", paste(missing_columns, collapse = ", ")))
}

duplicate_keys <- table$canonical_variant_key[duplicated(table$canonical_variant_key)]
validation_lines <- c(
  sprintf("reference_genome_build=%s", ifelse(is.null(reference_genome_build), "GRCh38", reference_genome_build)),
  sprintf("dry_run=%s", ifelse(is.null(dry_run), "true", dry_run)),
  sprintf("rows=%d", nrow(table)),
  sprintf("duplicate_keys=%d", length(unique(duplicate_keys))),
  sprintf("chromosomes=%s", paste(sort(unique(table$chromosome)), collapse = ","))
)

if (any(is.na(table$canonical_variant_key)) || any(table$canonical_variant_key == "")) {
  stop("Canonical variant keys must be present for every row.")
}

writeLines(validation_lines, con = output_path)