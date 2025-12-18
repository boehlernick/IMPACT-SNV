#!/usr/bin/env Rscript

# SeqArray Append Script for IMPACT-SNV
# This script adds additional annotations to GDS files:
# - Extracts and adds tier information from patho_score_calc
# - Creates boolean indicator nodes for ClinVar significance flags
# - Renames patho_score -> impact_score and patho_score_calc -> impact_score_calc
# - Optimizes the GDS file structure

# Usage: source("seqarray_append.R") and call process_gds_file(gds_path)
# Or run directly: Rscript seqarray_append.R <gds_file>

suppressPackageStartupMessages({
  library(SeqArray)
  library(gdsfmt)
  library(stringi)
})

#' Process a GDS file to add tier, clnsig flags, and rename score fields
#'
#' @param gds_path Path to the GDS file to process
#' @return NULL (modifies file in place)
process_gds_file <- function(gds_path) {
  if (!file.exists(gds_path)) {
    stop("GDS file not found: ", gds_path)
  }
  
  cat("Processing GDS file:", gds_path, "\n")
  
  # Open GDS file for writing
  gds <- seqOpen(gds_path, readonly = FALSE)
  
  # Ensure the file is closed when done or if an error occurs
  on.exit({
    if (exists("gds") && !is.null(gds)) {
      try(seqClose(gds), silent = TRUE)
    }
  }, add = TRUE)
  
  # Step 1: Add tier information extracted from patho_score_calc
  cat("Step 1: Adding tier information...\n")
  tryCatch({
    impact_score_calc <- seqGetData(gds, "annotation/info/patho_score_calc")
    tier <- as.integer(substr(impact_score_calc, 6, 6))
    seqAddValue(gds, "annotation/info/tier", tier, replace = TRUE)
    cat("  -> Added tier annotation (", sum(!is.na(tier) & tier > 0), " non-zero values)\n", sep = "")
  }, error = function(e) {
    warning("Could not add tier information: ", e$message)
  })
  
  # Step 2: Create ClinVar significance flag annotations
  cat("Step 2: Creating ClinVar significance boolean flags...\n")
  tryCatch({
    add_clnsig_flags(gds)
  }, error = function(e) {
    warning("Could not add clnsig flags: ", e$message)
  })
  
  # Step 3: Rename patho_score fields to impact_score
  cat("Step 3: Renaming score fields to impact_score...\n")
  tryCatch({
    rename_score_fields(gds)
  }, error = function(e) {
    warning("Could not rename score fields: ", e$message)
  })
  
  # Close before optimize
  seqClose(gds)
  gds <- NULL
  
  # Step 4: Optimize the GDS file
  cat("Step 4: Optimizing GDS file structure...\n")
  seqOptimize(gds_path)
  
  cat("Successfully processed:", gds_path, "\n\n")
}

#' Add ClinVar significance boolean flag annotations
#'
#' @param gds An open SeqArray GDS object (read/write mode)
add_clnsig_flags <- function(gds) {
  # Try to read clnsig data
  clnsig_raw <- tryCatch({
    seqGetData(gds, "annotation/info/FunctionalAnnotation/clnsig")
  }, error = function(e) {
    warning("No clnsig annotation found, skipping flag creation")
    return(NULL)
  })
  
  if (is.null(clnsig_raw)) return(invisible(NULL))
  
  clnsig_raw <- as.character(clnsig_raw)
  
  # Define the standard ClinVar flags to always include (label -> node)
  label_to_flag <- c(
    "Pathogenic"                           = "pathogenic",
    "Likely pathogenic"                    = "likely_pathogenic",
    "Uncertain significance"               = "uncertain_significance",
    "Likely benign"                        = "likely_benign",
    "Benign"                               = "benign",
    "Pathogenic, low penetrance"           = "pathogenic_low_penetrance",
    "Likely pathogenic, low penetrance"    = "likely_pathogenic_low_penetrance",
    "Established risk allele"              = "established_risk_allele",
    "Likely risk allele"                   = "likely_risk_allele",
    "Uncertain risk allele"                = "uncertain_risk_allele",
    "affects"                              = "affects",
    "association"                          = "association",
    "drug response"                        = "drug_response",
    "confers sensitivity"                  = "confers_sensitivity",
    "protective"                           = "protective",
    "other"                                = "other",
    "Conflicting Interpretations"          = "conflicting_interpretations_of_pathogenicity",
    "not provided"                         = "not_provided"
  )
  
  labels_lower <- tolower(names(label_to_flag))
  flags_nodes  <- unname(label_to_flag)
  
  # Be permissive for the "conflicting" label text
  labels_lower[labels_lower == "conflicting interpretations"] <-
    "conflicting interpretations of pathogenicity"
  
  cat("  Normalizing and tokenizing clnsig values...\n")
  
  # Identify missing/NA/empty entries
  trimmed    <- stri_trim_both(clnsig_raw)
  is_missing <- is.na(trimmed) | (trimmed == "") | (tolower(trimmed) == "na")
  
  # Build per-variant normalized tokens
  tokens_list <- lapply(seq_along(clnsig_raw), function(i) {
    x <- clnsig_raw[i]
    if (is_missing[i]) return(character(0))
    
    x <- tolower(x)
    x <- stri_replace_all_fixed(x, "_", " ")
    x <- stri_trim_both(x)
    
    # Protect phrase ", low penetrance"
    x <- stri_replace_all_regex(x, ",\\s*low\\s+penetrance", " ~~LP~~ ", vectorize_all = FALSE)
    
    # Unify delimiters to '|'
    x <- stri_replace_all_regex(x, "[/|;]", "|", vectorize_all = FALSE)
    x <- stri_replace_all_fixed(x, ",", "|",  vectorize_all = FALSE)
    x <- stri_replace_all_regex(x, "\\|+", "|", vectorize_all = FALSE)
    
    parts <- unlist(stri_split_fixed(x, "|", omit_empty = TRUE))
    parts <- stri_replace_all_fixed(parts, "~~LP~~", ", low penetrance")
    parts <- stri_trim_both(parts)
    parts <- stri_replace_all_regex(parts, "\\s+", " ")
    parts[nzchar(parts)]
  })
  
  cat("  Creating boolean indicator nodes for ClinVar flags...\n")
  
  # Get a handle to the "annotation/info" node
  info_node <- index.gdsn(gds, "annotation/info")
  
  # Ensure 'clnsig_flags' folder exists
  flags_node <- index.gdsn(info_node, "clnsig_flags", silent = TRUE)
  if (is.null(flags_node)) {
    flags_node <- addfolder.gdsn(info_node, "clnsig_flags")
    cat("  Created new folder: annotation/info/clnsig_flags\n")
  }
  
  # Build all flag vectors
  flag_vectors <- setNames(vector("list", length(flags_nodes)), flags_nodes)
  
  # Compute initial vectors from exact token membership
  for (j in seq_along(labels_lower)) {
    label_norm <- labels_lower[j]
    node_name  <- flags_nodes[j]
    
    bool_vector <- vapply(
      tokens_list,
      function(tokens) any(tokens == label_norm),
      logical(1L)
    )
    
    # 'not_provided' is TRUE for explicit token + also for missing/empty/"NA"
    if (node_name == "not_provided") {
      bool_vector <- bool_vector | is_missing
    }
    
    flag_vectors[[node_name]] <- bool_vector
  }
  
  # Enhance 'other': TRUE if explicit 'other' token OR any unrecognized token exists
  standard_set <- labels_lower
  has_tokens   <- lengths(tokens_list) > 0
  has_unrecognized <- vapply(tokens_list, function(tokens) any(!(tokens %in% standard_set)), logical(1L))
  
  if (!"other" %in% names(flag_vectors)) {
    flag_vectors[["other"]] <- has_unrecognized
  } else {
    flag_vectors[["other"]] <- flag_vectors[["other"]] | (has_tokens & has_unrecognized)
  }
  
  # Sanity check
  any_true <- Reduce(`|`, flag_vectors)
  if (any(!any_true)) {
    warning(sum(!any_true), " variants have no TRUE clnsig flag")
  }
  
  # Write all nodes to GDS (replace = TRUE for idempotency)
  for (node_name in names(flag_vectors)) {
    add.gdsn(flags_node, node_name, val = flag_vectors[[node_name]], storage = "bit1", replace = TRUE)
    cat(sprintf("    -> %-45s (%d hits)\n", node_name, sum(flag_vectors[[node_name]], na.rm = TRUE)))
  }
  
  cat("  ClinVar flag annotations complete.\n")
}

#' Rename patho_score fields to impact_score
#'
#' @param gds An open SeqArray GDS object (read/write mode)
rename_score_fields <- function(gds) {
  # Copy patho_score to impact_score, then delete patho_score
  tryCatch({
    patho_score <- seqGetData(gds, "annotation/info/patho_score")
    seqAddValue(gds, "annotation/info/impact_score", patho_score, replace = TRUE)
    seqDelete(gds, info.var = "patho_score")
    cat("  -> Renamed patho_score to impact_score\n")
  }, error = function(e) {
    cat("  -> patho_score not found or already renamed\n")
  })
  
  tryCatch({
    patho_score_calc <- seqGetData(gds, "annotation/info/patho_score_calc")
    seqAddValue(gds, "annotation/info/impact_score_calc", patho_score_calc, replace = TRUE)
    seqDelete(gds, info.var = "patho_score_calc")
    cat("  -> Renamed patho_score_calc to impact_score_calc\n")
  }, error = function(e) {
    cat("  -> patho_score_calc not found or already renamed\n")
  })
}

# Main execution when run as a script
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Usage: Rscript seqarray_append.R <gds_file> [<gds_file2> ...]\n")
    cat("  Processes one or more GDS files to add IMPACT annotations.\n")
    quit(status = 1)
  }
  
  for (gds_file in args) {
    tryCatch({
      process_gds_file(gds_file)
    }, error = function(e) {
      cat("Error processing", gds_file, ":", e$message, "\n")
    })
  }
}
