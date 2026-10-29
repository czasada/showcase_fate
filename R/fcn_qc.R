# Numbers behind the QC tab (section 2 of the WP1 degradomics report), aggregated from the final tables of all plates.
# QC rule of the run log: a run is flagged when its protein IDs fall below 85 % of its plate's median.

qc_min_share <- 0.85

#' QC metrics per run (protein IDs, summed raw intensity) with plate layout, run status and a running x position
#' over all plates (plate blocks of 96 in measurement order).
run_qc <- function(study) {
  study$intensities %>%
    dplyr::filter(!is.na(unnormalized_intensity)) %>%
    dplyr::summarise(
      n_ids = dplyr::n(), summed_intensity = sum(unnormalized_intensity),
      .by = c(plate_id, compound, replicate)
    ) %>%
    dplyr::left_join(
      study$plate_assignment %>% dplyr::select(plate_id, compound, replicate, sample, position, measurement_order),
      by = c("plate_id", "compound", "replicate")
    ) %>%
    dplyr::left_join(
      study$run_log %>% dplyr::select(plate_id, position, qc_status), by = c("plate_id", "position")
    ) %>%
    dplyr::mutate(
      share_of_median = n_ids / stats::median(n_ids),
      intensity_share_of_median = summed_intensity / stats::median(summed_intensity),
      .by = plate_id
    ) %>%
    dplyr::mutate(
      x = (match(plate_id, sort(unique(plate_id))) - 1) * 96 + measurement_order,
      well_type = dplyr::case_when(
        compound == vehicle_name ~ "DMSO",
        compound == positive_control ~ "positive control",
        .default = "compound"
      )
    ) %>%
    dplyr::arrange(x)
}

#' Box statistics of log2 intensities per run of one plate (precomputed, so the plot stays small).
#' @param type "unnormalized_intensity" or "normalized_intensity"
intensity_box_stats <- function(study, plate, type) {
  study$intensities %>%
    dplyr::filter(plate_id == plate, !is.na(.data[[type]])) %>%
    dplyr::mutate(value = log2(.data[[type]])) %>%
    dplyr::summarise(
      q1 = stats::quantile(value, 0.25), median = stats::median(value), q3 = stats::quantile(value, 0.75),
      lower = max(min(value), q1 - 1.5 * (q3 - q1)), upper = min(max(value), q3 + 1.5 * (q3 - q1)),
      .by = c(compound, replicate)
    ) %>%
    dplyr::left_join(
      study$plate_assignment %>%
        dplyr::filter(plate_id == plate) %>%
        dplyr::select(compound, replicate, sample, measurement_order),
      by = c("compound", "replicate")
    ) %>%
    dplyr::arrange(measurement_order)
}

#' Normalization factor per run: median log2 difference between raw and normalized intensity.
normalization_factors <- function(study) {
  factors <- study$intensities %>%
    dplyr::filter(!is.na(normalized_intensity)) %>%
    dplyr::summarise(
      factor = stats::median(log2(unnormalized_intensity) - log2(normalized_intensity)),
      .by = c(plate_id, compound, replicate)
    )
  run_qc(study) %>%
    dplyr::select(plate_id, compound, replicate, sample, x, well_type) %>%
    dplyr::left_join(factors, by = c("plate_id", "compound", "replicate"))
}

#' CV (%) of the DMSO wells per protein and plate, raw and normalized (proteins with at least min_values values).
dmso_cv <- function(study, min_values = 3) {
  study$intensities %>%
    dplyr::filter(compound == vehicle_name) %>%
    tidyr::pivot_longer(
      c(unnormalized_intensity, normalized_intensity),
      names_to = "type", values_to = "intensity"
    ) %>%
    dplyr::filter(!is.na(intensity)) %>%
    dplyr::summarise(
      cv = stats::sd(intensity) / mean(intensity) * 100, n = dplyr::n(),
      .by = c(plate_id, type, scout_id)
    ) %>%
    dplyr::filter(n >= min_values) %>%
    dplyr::mutate(type = ifelse(type == "normalized_intensity", "normalized", "raw"))
}

#' PCA of all runs on normalized log2 intensities of proteins measured in every run (centred, not scaled).
qc_pca <- function(study) {
  runs <- run_qc(study)
  m <- study$intensities %>%
    dplyr::mutate(run = paste(plate_id, compound, replicate, sep = "|"), value = log2(normalized_intensity)) %>%
    tidyr::pivot_wider(id_cols = scout_id, names_from = run, values_from = value) %>%
    tibble::column_to_rownames("scout_id") %>%
    as.matrix()
  m <- m[stats::complete.cases(m), , drop = FALSE]
  pca <- stats::prcomp(t(m), center = TRUE, scale. = FALSE)
  scores <- tibble::as_tibble(pca$x[, 1:3], rownames = "run") %>%
    tidyr::separate(run, c("plate_id", "compound", "replicate"), sep = "\\|") %>%
    dplyr::left_join(runs %>% dplyr::select(plate_id, compound, replicate, sample, well_type, qc_status),
      by = c("plate_id", "compound", "replicate")
    )
  list(scores = scores, variance = pca$sdev^2 / sum(pca$sdev^2) * 100, n_proteins = nrow(m))
}

#' Completeness per run (% of the plate's proteins measured) and valid values per protein and plate.
completeness <- function(study) {
  per_run <- run_qc(study) %>%
    dplyr::left_join(
      study$intensities %>% dplyr::summarise(n_proteins = dplyr::n_distinct(scout_id), .by = plate_id),
      by = "plate_id"
    ) %>%
    dplyr::mutate(completeness = n_ids / n_proteins * 100)
  per_protein <- study$intensities %>%
    dplyr::summarise(valid = sum(!is.na(unnormalized_intensity)), runs = dplyr::n(), .by = c(plate_id, scout_id)) %>%
    dplyr::count(plate_id, valid, runs)
  list(per_run = per_run, per_protein = per_protein)
}
