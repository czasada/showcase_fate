# Numbers behind the Proteins of interest tab: screen hits rolled up per protein, and one protein across all
# compounds. Hits use the Overview cutoffs (BH-adjusted p <= 0.05, log2 FC <= -1); the positive control is left out.

protein_categories <- c(established = "known neosubstrate", candidate = "candidate", none = "other")

#' One row per protein degraded by at least one compound: category, best reduction, number of compounds degrading
#' it and the best compound. Known neosubstrates first, then candidates, then by number of compounds.
protein_targets <- function(stats, identifier) {
  stats %>%
    dplyr::filter(hit, compound != positive_control) %>%
    dplyr::summarise(
      n_compounds = dplyr::n_distinct(compound),
      best_log2fc = min(log2_fold_change),
      best_compound = compound[which.min(log2_fold_change)],
      .by = c(scout_id, gene_name, neosubstrate_category)
    ) %>%
    dplyr::left_join(dplyr::select(identifier, scout_id, uniprot_id), by = "scout_id") %>%
    dplyr::mutate(
      category = factor(protein_categories[neosubstrate_category], levels = protein_categories),
      best_reduction = fold_change_to_degradation_efficiency(2^best_log2fc)
    ) %>%
    dplyr::arrange(category, dplyr::desc(n_compounds), dplyr::desc(best_reduction), gene_name)
}

#' Abundance rank per protein: median normalized log2 intensity over all DMSO wells, as "top x %" of the proteome.
abundance_rank <- function(study) {
  study$intensities %>%
    dplyr::filter(compound == vehicle_name, !is.na(normalized_intensity)) %>%
    dplyr::summarise(abundance = stats::median(log2(normalized_intensity)), .by = scout_id) %>%
    dplyr::mutate(top_percent = pmax(1, ceiling((1 - dplyr::percent_rank(abundance)) * 100)))
}

#' Minimum detectable change of one protein, as log2 fold change: -2 x SD of its log2 intensity over the DMSO wells
#' (SD per plate, averaged over plates). Drops smaller than this are within the protein's DMSO noise.
min_detectable_change <- function(study, protein_id) {
  sds <- study$intensities %>%
    dplyr::filter(scout_id == protein_id, compound == vehicle_name, !is.na(normalized_intensity)) %>%
    dplyr::summarise(sd = stats::sd(log2(normalized_intensity)), .by = plate_id)
  -2 * mean(sds$sd, na.rm = TRUE)
}

#' One protein across all compounds: log2 FC, p-value, hit call, plate and how many proteins each compound degrades.
protein_profile <- function(stats, protein_id) {
  n_degraded <- compound_summary(stats) %>% dplyr::select(compound, n_degraded)
  stats %>%
    dplyr::filter(scout_id == protein_id, compound != positive_control, !is.na(log2_fold_change)) %>%
    dplyr::left_join(n_degraded, by = "compound") %>%
    dplyr::arrange(compound) %>%
    dplyr::mutate(order = dplyr::row_number())
}

# Value shown for a log2 fold change under the metric switch.
protein_metric_value <- function(log2fc, metric) {
  switch(metric,
    log2fc = log2fc,
    reduction = fold_change_to_degradation_efficiency(2^log2fc),
    remaining = 2^log2fc * 100
  )
}
protein_metric_label <- c(log2fc = "Log2 fold change", reduction = "Reduction (%)", remaining = "Remaining (%)")
