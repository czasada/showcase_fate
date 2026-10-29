# Numbers behind the Statistics tab (section 3 of the WP1 degradomics report): each compound vs the DMSO wells of its
# plate, as delivered in the final tables. Default cutoffs match the Overview (BH-adjusted p <= 0.05, log2 FC <= -1).

# Volcano axes as in the 0609 volcano report: x from -6 to 6 when all proteins fit, else -10 to 10 (proteins beyond
# are drawn on the edge); y limits extended to at least these values so volcanos stay comparable.
volcano_limits <- list(min_p = 1e-10, min_standard_error = 10^-1.5, max_standard_error = 10^0.5)
volcano_x_limit <- function(log2_fold_changes) if (all(abs(log2_fold_changes) <= 6, na.rm = TRUE)) 6 else 10

#' Statistics with regulation call under the given cutoffs and the standard error of each log2 fold change
#' (|log2 fold change / t-statistic|, not defined where t is 0 or missing; as in the 0609 volcano report).
#' @param stats output of annotated_statistics()
regulated_statistics <- function(stats, p_col = "p_value_corrected", p_cutoff = hit_p_cutoff,
                                 log2fc_cutoff = abs(hit_log2fc_cutoff), direction = "only_down") {
  stats %>%
    dplyr::filter(!is.na(log2_fold_change)) %>%
    add_regulation(p_col, p_cutoff, log2fc_cutoff, direction) %>%
    dplyr::mutate(
      standard_error = dplyr::if_else(!is.na(t_statistic) & t_statistic != 0, abs(log2_fold_change / t_statistic), NA),
      has_imputed = grepl("i", valid_value_string)
    )
}

#' Hits per compound: down / up, established neosubstrates among the down hits, tested proteins and the share of
#' hits involving imputed values (report: data completeness and imputed values).
compound_hit_counts <- function(reg) {
  reg %>%
    dplyr::summarise(
      tested = dplyr::n(),
      down = sum(regulation == "down"),
      up = sum(regulation == "up"),
      neo_down = sum(regulation == "down" & neosubstrate_category == "established"),
      imputed_share = ifelse(sum(regulation != "none") > 0, mean(has_imputed[regulation != "none"]) * 100, NA),
      .by = c(compound, plate_id)
    ) %>%
    dplyr::arrange(dplyr::desc(down + up), compound)
}

#' Protein shown next to the volcano before anything is clicked: the compound's strongest hit (lowest p-value among
#' the regulated proteins), else its most reduced protein.
default_protein <- function(reg_compound, p_col = "p_value_corrected") {
  hits <- reg_compound %>% dplyr::filter(regulation != "none")
  pick <- if (nrow(hits) > 0) {
    dplyr::slice_min(hits, .data[[p_col]], n = 1, with_ties = FALSE)
  } else {
    dplyr::slice_min(reg_compound, log2_fold_change, n = 1, with_ties = FALSE)
  }
  pick$scout_id
}

#' The n other compounds that degrade a protein most strongly under the current cutoffs.
top_degraders <- function(reg, protein_id, exclude, n = 3) {
  reg %>%
    dplyr::filter(scout_id == protein_id, regulation == "down", compound != exclude) %>%
    dplyr::slice_min(log2_fold_change, n = n, with_ties = FALSE) %>%
    dplyr::pull(compound)
}

#' Replicate values of one protein for the given compounds, grouped by plate: the plate's DMSO wells first, then its
#' compounds in the given order (plates in order of their first compound). Compound values are the imputed
#' intensities the statistics used (imputed ones flagged); DMSO is not imputed, so missing DMSO values are left out.
#' @return tibble(group, role = DMSO / compound / other compound, plate_id, replicate, value, is_imputed)
protein_replicates <- function(study, protein_id, compounds) {
  plates <- study$statistics %>%
    dplyr::distinct(compound, plate_id) %>%
    dplyr::filter(compound %in% compounds) %>%
    dplyr::arrange(match(compound, compounds))
  plate_order <- unique(plates$plate_id)
  study$intensities %>%
    dplyr::filter(scout_id == protein_id, plate_id %in% plate_order) %>%
    dplyr::left_join(dplyr::rename(plates, compound_plate = plate_id), by = "compound") %>%
    dplyr::filter(compound == vehicle_name | compound_plate == plate_id) %>%
    dplyr::mutate(
      is_imputed = is.na(normalized_intensity),
      value = log2(ifelse(compound == vehicle_name, normalized_intensity, imputed_intensity)),
      role = dplyr::case_when(
        compound == vehicle_name ~ "DMSO", compound == compounds[1] ~ "compound", .default = "other compound"
      ),
      group = ifelse(role == "DMSO", paste0("DMSO (plate ", sub("P", "", plate_id), ")"), compound)
    ) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::arrange(match(plate_id, plate_order), role != "DMSO", match(compound, compounds), replicate) %>%
    dplyr::select(group, role, plate_id, replicate, value, is_imputed)
}
