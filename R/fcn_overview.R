# Numbers behind the Overview tab, aggregated from the final tables and the run log of read_study().
# Wording targets readers outside proteomics (reduction in %, plain stage names).

# Hit = degraded: BH-adjusted p <= 0.05 and log2 FC <= -1 (at least 50 % reduction). Demo convention.
hit_p_cutoff <- 0.05
hit_log2fc_cutoff <- -1
pos_ctrl_acceptance <- 75

#' Statistics with gene names, annotation and hit call.
annotated_statistics <- function(study, p_cutoff = hit_p_cutoff, log2fc_cutoff = hit_log2fc_cutoff) {
  genes <- study$identifier %>% dplyr::select(scout_id, gene_name)
  study$statistics %>%
    dplyr::left_join(genes, by = "scout_id") %>%
    dplyr::left_join(study$annotations, by = "scout_id") %>%
    dplyr::mutate(
      hit = !is.na(p_value_corrected) & p_value_corrected <= p_cutoff & log2_fold_change <= log2fc_cutoff
    )
}

#' One row per tested compound (positive control excluded): degraded proteins, neosubstrate and candidate hits.
compound_summary <- function(stats) {
  stats %>%
    dplyr::filter(compound != positive_control) %>%
    dplyr::summarise(
      n_degraded = sum(hit),
      n_neo = sum(hit & neosubstrate_category == "established"),
      n_candidates = sum(hit & neosubstrate_category == "candidate"),
      .by = c(compound, plate_id)
    )
}

#' Positive control readout (GSPT1 reduction, %) per plate from the statistics.
positive_control_per_plate <- function(stats) {
  stats %>%
    dplyr::filter(compound == positive_control, gene_name == positive_control_readout) %>%
    dplyr::transmute(plate_id, reduction = fold_change_to_degradation_efficiency(2^log2_fold_change))
}

#' Positive control readout per replicate: GSPT1 intensity relative to the mean of the plate's DMSO wells.
positive_control_replicates <- function(study) {
  readout_id <- study$identifier$scout_id[study$identifier$gene_name == positive_control_readout]
  study$intensities %>%
    dplyr::filter(scout_id %in% readout_id, compound %in% c(vehicle_name, positive_control)) %>%
    dplyr::mutate(vehicle_mean = mean(normalized_intensity[compound == vehicle_name], na.rm = TRUE), .by = plate_id) %>%
    dplyr::filter(compound == positive_control) %>%
    dplyr::transmute(
      plate_id, replicate,
      reduction = fold_change_to_degradation_efficiency(normalized_intensity / vehicle_mean)
    )
}

plate_status <- function(study) {
  study$run_log %>%
    dplyr::summarise(
      acquired = sum(!is.na(acquired_date)), n_runs = dplyr::n(), week = min(planned_date),
      .by = plate_id
    ) %>%
    dplyr::mutate(
      status = dplyr::case_when(
        plate_id %in% study$statistics$plate_id ~ "results available",
        acquired > 0 ~ "acquiring",
        .default = "planned"
      ),
      week = paste("wk", as.integer(week - min(week)) %/% 7 + 1)
    )
}

#' Tiles of the progress section.
progress_summary <- function(study, stats) {
  plates <- plate_status(study)
  log <- study$run_log
  n_passed <- sum(log$qc_status == "passed")
  n_flagged <- sum(log$qc_status == "flagged")
  compounds <- compound_summary(stats)
  planned <- setdiff(unique(log$compound), c(vehicle_name, positive_control))
  list(
    n_plates = nrow(plates),
    n_plates_done = sum(plates$status == "results available"),
    runs_acquired = sum(!is.na(log$acquired_date)),
    runs_total = nrow(log),
    qc_pass_rate = n_passed / (n_passed + n_flagged) * 100,
    n_qc_flagged = n_flagged,
    pos_ctrl = stats::median(positive_control_per_plate(stats)$reduction),
    n_compounds_planned = length(planned),
    n_compounds_results = nrow(compounds),
    n_neo_compounds = sum(compounds$n_neo > 0)
  )
}

#' Cumulative runs per day: plan (all runs) and acquired (so far).
runs_against_plan <- function(run_log) {
  cumulative <- function(dates) {
    tibble::tibble(date = sort(dates)) %>%
      dplyr::count(date) %>%
      dplyr::mutate(runs = cumsum(n))
  }
  list(
    plan = cumulative(run_log$planned_date),
    acquired = cumulative(run_log$acquired_date[!is.na(run_log$acquired_date)])
  )
}

#' One row per plate for the plates table.
plates_table <- function(study, stats) {
  dmso_cv <- study$intensities %>%
    dplyr::filter(compound == vehicle_name, !is.na(normalized_intensity)) %>%
    dplyr::summarise(
      cv = stats::sd(normalized_intensity) / mean(normalized_intensity) * 100, n = dplyr::n(),
      .by = c(plate_id, scout_id)
    ) %>%
    dplyr::filter(n >= 3) %>%
    dplyr::summarise(dmso_cv = stats::median(cv), .by = plate_id)
  qc <- study$run_log %>%
    dplyr::summarise(
      runs_passed = sum(qc_status == "passed"), runs_checked = sum(qc_status %in% c("passed", "flagged")),
      .by = plate_id
    )
  neo <- compound_summary(stats) %>%
    dplyr::summarise(n_neo_compounds = sum(n_neo > 0), .by = plate_id)
  plate_status(study) %>%
    dplyr::left_join(positive_control_per_plate(stats), by = "plate_id") %>%
    dplyr::left_join(dmso_cv, by = "plate_id") %>%
    dplyr::left_join(qc, by = "plate_id") %>%
    dplyr::left_join(neo, by = "plate_id")
}

#' Thresholds of the card badges. These are demo conventions, not literature values.
#'
#' @param broad_min known neosubstrates degraded for "broad" (otherwise the degraded ones are named)
#' @param selective_share share of neosubstrates among all degraded proteins for "selective"
#' @param unselective_share share below which a compound with at least unselective_min_degraded degraded
#'   proteins gets the caution badge
card_badge_rules <- function(broad_min = 3, selective_share = 0.25, unselective_share = 0.1,
                             unselective_min_degraded = 10) {
  list(
    broad_min = broad_min, selective_share = selective_share, unselective_share = unselective_share,
    unselective_min_degraded = unselective_min_degraded
  )
}

#' Profile badge: which neosubstrates does the compound degrade?
#' "Broad" counts established neosubstrates only; below that, all degraded neosubstrates and candidates are named.
#' @param neo_hits gene names of the degraded established neosubstrates and candidates
#' @param zinc_finger logical per neo_hits
#' @param established logical per neo_hits (FALSE = candidate)
profile_badge <- function(neo_hits, zinc_finger, n_degraded, established = rep(TRUE, length(neo_hits)),
                          rules = card_badge_rules()) {
  n_neo <- length(neo_hits)
  if (sum(established) >= rules$broad_min) {
    label <- if (mean(zinc_finger) >= 0.5) "Broad zinc-finger degrader" else "Broad neosubstrate degrader"
    list(label = label, level = "strong")
  } else if (n_neo > 0) {
    list(label = paste(paste(neo_hits, collapse = " / "), "degrader"), level = "moderate")
  } else if (n_degraded > 0) {
    list(label = "No neosubstrate degraded", level = "neutral")
  } else {
    list(label = "Weak effect", level = "neutral")
  }
}

#' Selectivity badge from the share of neosubstrates among all degraded proteins (NULL when in between).
selectivity_badge <- function(n_neo, n_degraded, rules = card_badge_rules()) {
  if (n_degraded == 0 || n_neo == 0) {
    return(NULL)
  }
  share <- n_neo / n_degraded
  if (n_degraded >= rules$unselective_min_degraded && share < rules$unselective_share) {
    list(label = "Many other proteins reduced", level = "caution")
  } else if (share >= rules$selective_share) {
    list(label = "Selective", level = "ok")
  }
}

#' Compounds on the card wheel: all that degrade at least one established neosubstrate, most neosubstrates first,
#' then fewest degraded proteins overall.
card_compounds <- function(stats) {
  compound_summary(stats) %>%
    dplyr::filter(n_neo > 0) %>%
    dplyr::arrange(dplyr::desc(n_neo + n_candidates), n_degraded, compound) %>%
    dplyr::pull(compound)
}

#' Data of one compound card. The reduction list shows the strongest-reduced neosubstrates and candidates,
#' including non-significant ones (hit = FALSE).
#' @return list(compound, plate_id, wells, n_degraded, n_neo, n_candidates, badges, top)
compound_card <- function(cmpd, stats, plate_assignment, n_top = 4, rules = card_badge_rules()) {
  cmpd_stats <- stats %>% dplyr::filter(compound == cmpd)
  neo_stats <- cmpd_stats %>% dplyr::filter(neosubstrate_category != "none")
  hits <- neo_stats %>% dplyr::filter(hit)
  n_degraded <- sum(cmpd_stats$hit)
  top <- neo_stats %>%
    dplyr::filter(log2_fold_change < 0) %>%
    dplyr::slice_min(log2_fold_change, n = n_top, with_ties = FALSE) %>%
    dplyr::transmute(gene_name, hit, reduction = round(fold_change_to_degradation_efficiency(2^log2_fold_change)))
  wells <- plate_assignment %>% dplyr::filter(compound == cmpd)
  badges <- list(
    profile_badge(hits$gene_name, hits$zinc_finger, n_degraded, hits$neosubstrate_category == "established", rules),
    selectivity_badge(nrow(hits), n_degraded, rules)
  )
  list(
    compound = cmpd, plate_id = cmpd_stats$plate_id[1], wells = paste(wells$position, collapse = " / "),
    n_degraded = n_degraded, n_neo = sum(hits$neosubstrate_category == "established"),
    n_candidates = sum(hits$neosubstrate_category == "candidate"), badges = purrr::compact(badges), top = top
  )
}
