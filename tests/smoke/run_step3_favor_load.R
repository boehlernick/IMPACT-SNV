
# tests/smoke/run_step3_favor_load.R
box::use(
  modules/steps/step3_favor_load[run],
  readr[write_csv],
  dplyr[glimpse],
  fs
)

main <- function() {
  fixture_dir <- "tests/smoke/fixtures"
  out_path <- fs::path(fixture_dir, "favor_step3_loaded.csv")

  df <- run(fixture_dir = fixture_dir)
  dplyr::glimpse(df)

  readr::write_csv(df, out_path)
  message(sprintf("Wrote combined FAVOR to: %s (n=%d)", out_path, nrow(df)))
}

main()
