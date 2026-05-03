#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeqArray)
  library(gdsfmt)
})

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

normalize_chromosome <- function(chromosome) {
  chromosome <- toupper(trimws(as.character(chromosome)))
  chromosome <- sub("^CHR", "", chromosome)
  if (chromosome %in% c(as.character(1:22), "X", "Y")) {
    return(chromosome)
  }
  if (chromosome %in% c("23", "24")) {
    stop("Numeric sex chromosome aliases 23 and 24 are not supported in this scaffold.")
  }
  stop(sprintf("Unsupported chromosome value: %s", chromosome))
}

extract_allele_pair <- function(allele_value) {
  if (is.matrix(allele_value) || is.data.frame(allele_value)) {
    allele_value <- apply(as.matrix(allele_value), 1, function(row) paste(row[nzchar(row)], collapse = ","))
  }
  allele_value <- as.character(allele_value)
  allele_value <- gsub("\\s+", "", allele_value)

  parse_one <- function(value) {
    pieces <- unlist(strsplit(value, "[,/|]"))
    pieces <- pieces[nzchar(pieces)]
    if (length(pieces) == 0) {
      return(c(NA_character_, NA_character_))
    }
    ref <- pieces[[1]]
    alt <- if (length(pieces) > 1) paste(pieces[-1], collapse = ",") else NA_character_
    c(ref, alt)
  }

  t(vapply(allele_value, parse_one, character(2)))
}

main <- function() {
  args <- parse_args()
  input_path <- args[["input"]]
  input_type <- args[["input-type"]]
  output_path <- args[["output"]]
  reference_genome_build <- args[["reference-genome-build"]]

  if (is.null(input_path) || is.null(input_type) || is.null(output_path)) {
    stop("Usage: --input <path> --input-type <gds|vcf> --output <path> [--reference-genome-build <build>]")
  }

  if (input_type == "vcf") {
    con <- if (grepl("\\.gz$", input_path, ignore.case = TRUE)) gzfile(input_path, open = "rt") else file(input_path, open = "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)

    variant_rows <- list()
    variant_counter <- 0L
    while (TRUE) {
      line <- readLines(con, n = 1)
      if (length(line) == 0) {
        break
      }
      if (startsWith(line, "#")) {
        next
      }

      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) < 5) {
        next
      }

      chromosome <- normalize_chromosome(fields[[1]])
      position <- fields[[2]]
      ref <- fields[[4]]
      alt_values <- unlist(strsplit(fields[[5]], ",", fixed = TRUE))
      alt_values <- alt_values[nzchar(alt_values)]

      for (alt in alt_values) {
        variant_counter <- variant_counter + 1L
        variant_rows[[variant_counter]] <- data.frame(
          variant.id = variant_counter,
          chromosome = chromosome,
          position = as.integer(position),
          ref = ref,
          alt = alt,
          canonical_variant_key = paste(chromosome, position, ref, alt, sep = "-"),
          reference_genome_build = reference_genome_build,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(variant_rows) == 0) {
      variant_table <- data.frame(
        variant.id = integer(),
        chromosome = character(),
        position = integer(),
        ref = character(),
        alt = character(),
        canonical_variant_key = character(),
        reference_genome_build = character(),
        stringsAsFactors = FALSE
      )
    } else {
      variant_table <- do.call(rbind, variant_rows)
    }
  } else {
    gds <- seqOpen(input_path, readonly = TRUE)
    on.exit(try(seqClose(gds), silent = TRUE), add = TRUE)

    variant_id <- seqGetData(gds, "variant.id")
    chromosome <- seqGetData(gds, "chromosome")
    position <- seqGetData(gds, "position")

    allele <- tryCatch(seqGetData(gds, "allele"), error = function(e) NULL)
    if (is.null(allele)) {
      stop("Unable to read allele information from GDS; canonical identity extraction requires allele data.")
    }

    allele_pairs <- extract_allele_pair(allele)
    ref <- allele_pairs[, 1]
    alt <- allele_pairs[, 2]

    normalized_chromosome <- vapply(chromosome, normalize_chromosome, character(1))
    canonical_variant_key <- paste(normalized_chromosome, position, ref, alt, sep = "-")

    variant_table <- data.frame(
      variant.id = variant_id,
      chromosome = normalized_chromosome,
      position = position,
      ref = ref,
      alt = alt,
      canonical_variant_key = canonical_variant_key,
      reference_genome_build = reference_genome_build,
      stringsAsFactors = FALSE
    )
  }

  write.table(variant_table, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
}

main()