# Summaries and hit calls on the synthetic data, following the WP1 degradomics report.

#' Regulation call per protein and comparison: "down", "up" or "none".
#' Cutoffs are inclusive (p <= cutoff, |log2 FC| >= cutoff).
add_regulation <- function(statistics, p_col = "p_value_corrected", p_cutoff = 0.05, log2fc_cutoff = 1,
                           direction = "both-sided") {
  statistics %>%
    dplyr::mutate(regulation = dplyr::case_when(
      .data[[p_col]] > p_cutoff ~ "none",
      direction != "only_down" & log2_fold_change >= log2fc_cutoff ~ "up",
      direction != "only_up" & log2_fold_change <= -log2fc_cutoff ~ "down",
      .default = "none"
    ))
}

# same name as in the report templates
fold_change_to_degradation_efficiency <- function(fc, as_percent = TRUE, clamp = TRUE) { # nolint: object_length_linter.
  eff <- 1 - fc
  if (clamp) eff <- pmax(eff, 0) # fc > 1 (enrichment) -> 0%
  if (as_percent) eff * 100 else eff
}
