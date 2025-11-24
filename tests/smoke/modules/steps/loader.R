
# tests/smoke/modules/steps/loader.R

box::use(
  readr[read_tsv, read_csv, cols, col_character],
  dplyr[bind_rows, select, mutate, across],
  tibble[tibble],
  fs,
  purrr[map]
)

# ---- helpers: suffix-based typing, aligned with Step3-favorannotator-rap/favorannotator.R ----

to_logical <- function(x) {
  x <- trimws(tolower(x))
  dplyr::case_when(
    x %in% c("true", "t", "1", "yes", "y")  ~ TRUE,
    x %in% c("false", "f", "0", "no", "n")  ~ FALSE,
    TRUE                                    ~ NA
  )
}

coerce_by_suffix_one <- function(df) {
  nms <- names(df)
  str_cols <- grep("\\.String$",  nms, value = TRUE)
  val_cols <- grep("\\.Valid$",   nms, value = TRUE)
  int_cols <- grep("\\.Int32$",   nms, value = TRUE)
  dbl_cols <- grep("\\.Float64$", nms, value = TRUE)
  
  if (length(str_cols)) {
    df <- dplyr::mutate(df, dplyr::across(dplyr::all_of(str_cols), as.character))
  }
  if (length(val_cols)) {
    df <- dplyr::mutate(df, dplyr::across(dplyr::all_of(val_cols), to_logical))
  }
  if (length(int_cols)) {
    suppressWarnings({
      df <- dplyr::mutate(df, dplyr::across(dplyr::all_of(int_cols), ~ as.integer(.x)))
    })
  }
  if (length(dbl_cols)) {
    suppressWarnings({
      df <- dplyr::mutate(df, dplyr::across(dplyr::all_of(dbl_cols), ~ as.numeric(.x)))
    })
  }
  df
}

# ---- main API ----

load_favor_subsets <- function(map_path,
                               required_cols = c(
                                 "Chromosome.String", "Position.String", "RefVcf.String", "AltVcf.String",
                                 "VariantVcf.String", "MetasvmPred.String", "CagePromoter.String",
                                 "Genehancer.String", "Linsight.Float64", "CaddPhred.Float64"
                               )) {
  if (!fs::file_exists(map_path)) {
    stop(sprintf("FAVOR map not found: %s", map_path))
  }
  
  map_tbl <- readr::read_tsv(map_path, col_names = c("chrom", "csv_path"), show_col_types = FALSE)
  map_tbl <- dplyr::mutate(map_tbl, exists = fs::file_exists(csv_path))
  
  missing_files <- dplyr::filter(map_tbl, !exists)
  if (nrow(missing_files) > 0) {
    warning(sprintf("Skipping %d missing FAVOR CSV(s):\n- %s",
                    nrow(missing_files), paste(missing_files$csv_path, collapse = "\n- ")))
  }
  present <- dplyr::filter(map_tbl, exists)
  
  # --- KEY: read every shard as character, then coerce by suffix *before* binding ---
  load_one <- function(path) {
    df <- readr::read_csv(path,
                          show_col_types = FALSE,
                          col_types = readr::cols(.default = readr::col_character()))
    df <- coerce_by_suffix_one(df)
    
    # ensure keys are character (they are *.String)
    key_cols <- c("Chromosome.String", "Position.String", "RefVcf.String", "AltVcf.String")
    key_cols <- intersect(key_cols, names(df))
    if (length(key_cols)) {
      df <- dplyr::mutate(df, dplyr::across(dplyr::all_of(key_cols), as.character))
    }
    df
  }
  
  pieces <- purrr::map(present$csv_path, load_one)
  
  if (length(pieces) == 0) {
    return(tibble())
  }
  
  # Optional debug: verify classes are uniform (uncomment if needed)
  # \print(lapply(pieces, function(x) sapply(x[c("RefVcf.String","AltVcf.String")], class)))
  
  
  all <- dplyr::bind_rows(pieces, .id = ".file_index")
  
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, names(all))
    if (length(missing_cols) > 0) {
      warning(sprintf(
        "FAVOR subset is missing %d expected columns:\n- %s",
        length(missing_cols),
        paste(missing_cols, collapse = "\n- ")
      ))
    }
  }
  
  all
}

# export for box
box::export(load_favor_subsets)
