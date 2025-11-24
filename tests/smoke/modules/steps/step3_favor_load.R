
# modules/steps/step3_favor_load.R
box::use(
  ./loader[load_favor_subsets],
  tibble[tibble],
  dplyr[mutate],
  fs
)

#' Step 3: Load FAVOR locally from subset map
#'
#' @param fixture_dir Directory holding `favor_subset_map.txt` produced by build_favor_subset.py
#' @param map_name Optional name of map file (default: "favor_subset_map.txt")
#' @return A tibble with combined FAVOR annotations for smoke test variants.
run <- function(fixture_dir, map_name = "favor_subset_map.txt") {
  map_path <- fs::path(fixture_dir, map_name)
  df <- load_favor_subsets(map_path = map_path)

  # Optionally decorate with provenance for debugging
  df <- dplyr::mutate(df, .source = map_path)

  df
}
