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

#' Per-sample fold change to the mean vehicle intensity (imputed values included and flagged).
sample_fold_changes <- function(intensities) {
  intensities %>%
    dplyr::mutate(intensity = 2^imputed_log2_intensity) %>%
    dplyr::mutate(mean_vehicle = mean(intensity[compound == vehicle_name]), .by = protein_id) %>%
    dplyr::mutate(
      fold_change = intensity / mean_vehicle,
      log2_fold_change = log2(fold_change),
      degradation_efficiency_percent = fold_change_to_degradation_efficiency(fold_change)
    ) %>%
    dplyr::select(protein_id, sample, compound, replicate, log2_fold_change, degradation_efficiency_percent, is_imputed)
}

#' ID count, summed raw intensity and imputed count per sample.
sample_qc_summary <- function(intensities, samples) {
  intensities %>%
    dplyr::summarise(
      n_ids = sum(!is.na(raw_log2_intensity)),
      summed_intensity = sum(2^raw_log2_intensity, na.rm = TRUE),
      n_imputed = sum(is_imputed),
      .by = sample
    ) %>%
    dplyr::left_join(samples, by = "sample")
}

#' Coefficient of variation (%) per protein and condition on measured, normalized intensities.
cv_per_condition <- function(intensities) {
  intensities %>%
    dplyr::filter(!is_imputed) %>%
    dplyr::mutate(intensity = 2^normalized_log2_intensity) %>%
    dplyr::summarise(
      cv = stats::sd(intensity) / mean(intensity) * 100,
      n = dplyr::n(),
      .by = c(protein_id, compound)
    ) %>%
    dplyr::filter(n >= 2)
}

#' PCA on proteins without missing values (normalized log2 intensities, centred, not scaled).
pca_scores <- function(intensities, samples) {
  m <- intensities %>%
    tidyr::pivot_wider(id_cols = protein_id, names_from = sample, values_from = normalized_log2_intensity) %>%
    tibble::column_to_rownames("protein_id") %>%
    as.matrix()
  m <- m[stats::complete.cases(m), , drop = FALSE]
  pca <- stats::prcomp(t(m), center = TRUE, scale. = FALSE)
  variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
  scores <- tibble::as_tibble(pca$x[, 1:3], rownames = "sample") %>%
    dplyr::left_join(samples, by = "sample")
  list(scores = scores, variance = variance, n_proteins = nrow(m))
}

#' Pairwise overlap of regulated proteins between comparisons (same direction).
pairwise_overlap <- function(regulated, compounds) {
  sets <- split(paste(regulated$protein_id, regulated$regulation), factor(regulated$compound, levels = compounds))
  tidyr::expand_grid(compound_a = compounds, compound_b = compounds) %>%
    dplyr::mutate(n_shared = purrr::map2_int(compound_a, compound_b, \(a, b) length(intersect(sets[[a]], sets[[b]]))))
}

#' Intersections of regulated proteins: which combination of compounds each protein is a hit in.
hit_intersections <- function(regulated, proteins) {
  regulated %>%
    dplyr::summarise(
      compounds = paste(sort(unique(compound)), collapse = " & "),
      n_compounds = dplyr::n_distinct(compound),
      .by = c(protein_id, regulation)
    ) %>%
    dplyr::left_join(proteins %>% dplyr::select(protein_id, gene_name, neosubstrate_category), by = "protein_id")
}
