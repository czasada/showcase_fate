# Synthetic degradomics screen for the public demo dashboard.
# Everything here is generated from a fixed seed - no project data is read. The output mimics what the pipeline
# delivers: one flat-header final table per plate (read with prepare_final_table()), a sample annotation per plate
# (read with infer_plate_layout()) and a run log with plan, acquisition dates and QC status.

# Established molecular-glue neosubstrates from the public literature. Background proteins are
# synthetic (SYNG....), candidates are synthetic zinc-finger-like entries (SYNZF..).
established_neosubstrates <- c(
  "IKZF1", "IKZF2", "IKZF3", "GSPT1", "GSPT2", "CSNK1A1", "SALL4", "ZFP91",
  "ZNF692", "ZNF827", "WIZ", "FIZ1", "ZBTB16", "E4F1", "RAB28", "PDE6D"
)
zinc_finger_neosubstrates <- c(
  "IKZF1", "IKZF2", "IKZF3", "SALL4", "ZFP91", "ZNF692", "ZNF827", "WIZ", "FIZ1", "ZBTB16", "E4F1"
)

vehicle_name <- "DMSO"
positive_control <- "POS-CTRL"
positive_control_readout <- "GSPT1"
n_plates <- 3
n_measured_plates <- 2
compounds_per_plate <- 44

# Planted compounds with a known profile (log2 FC per neosubstrate, number of other proteins down).
compound_profiles <- list(
  "CMPD-0007" = list(neo = c(IKZF1 = -2.6, IKZF3 = -2.9, ZFP91 = -1.4, ZNF692 = -1.1), n_off = 8),
  "CMPD-0019" = list(neo = c(GSPT1 = -2.3, GSPT2 = -1.5), n_off = 30),
  "CMPD-0034" = list(neo = c(CSNK1A1 = -2.1, IKZF1 = -1.6, IKZF3 = -1.2), n_off = 6),
  "CMPD-0052" = list(
    neo = c(
      IKZF1 = -2.2, IKZF2 = -1.9, IKZF3 = -2.4, ZFP91 = -1.7, ZNF692 = -1.5, SALL4 = -1.3,
      WIZ = -1.8, FIZ1 = -1.2, SYNZF03 = -1.6, SYNZF07 = -1.3
    ),
    n_off = 12
  ),
  "CMPD-0071" = list(neo = c(PDE6D = -0.8), n_off = 2),
  "POS-CTRL" = list(neo = c(GSPT1 = -2.6, GSPT2 = -1.2), n_off = 10)
)

# How often each neosubstrate comes up in the drawn compounds (relative weights, made up but ordered like the
# literature: GSPT1 and the Ikaros family are the most frequent CRBN neosubstrates).
neosubstrate_weights <- c(
  GSPT1 = 8, IKZF1 = 7, IKZF3 = 6, CSNK1A1 = 4, ZFP91 = 3, IKZF2 = 3, GSPT2 = 2, SALL4 = 2, ZNF692 = 2,
  WIZ = 1.5, FIZ1 = 1, ZNF827 = 1, ZBTB16 = 0.7, E4F1 = 0.5, RAB28 = 0.5, PDE6D = 0.5
)
candidate_weights <- c(SYNZF03 = 3, SYNZF07 = 2, SYNZF01 = 1, SYNZF05 = 1)

#' Write the synthetic study to a folder: final_table_<plate>_flat_header.csv and sample_annotation_<plate>.csv for
#' each measured plate, plus run_log.csv for all plates.
#' @return the folder
write_synthetic_study <- function(folder, seed = 640, n_proteins = 4000, today = Sys.Date()) {
  set.seed(seed)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  proteins <- make_proteins(n_proteins)
  layout <- make_plate_layout()
  run_log <- make_run_log(layout, today)
  # three compound runs are under-loaded; the QC rule (protein IDs < 85 % of the plate median) flags them. Flagged
  # runs stay in the final table, as in the pipeline.
  can_fail <- !is.na(run_log$acquired_date) &
    !run_log$compound %in% c(vehicle_name, positive_control, names(compound_profiles))
  underloaded <- run_log$sample_id[sample(which(can_fail), 3)]

  qc <- list()
  for (i in seq_len(n_measured_plates)) {
    plate <- unique(layout$plate_id)[i]
    plate_layout <- layout %>% dplyr::filter(plate_id == plate)
    # small plate effect, cancelled out by comparing each compound to the DMSO wells of its own plate
    plate_data <- make_final_table(proteins, plate_layout, underloaded, plate_offset = 0.25 * (i - 1))
    qc[[plate]] <- plate_data$qc
    data.table::fwrite(plate_data$final_table, file.path(folder, paste0("final_table_", plate, "_flat_header.csv")))
    data.table::fwrite(
      plate_layout %>% dplyr::select(experiment_id, compound, replicate, cell_line, compound_concentration, time),
      file.path(folder, paste0("sample_annotation_", plate, ".csv"))
    )
  }
  run_log <- run_log %>%
    dplyr::left_join(dplyr::bind_rows(qc), by = "sample_id") %>%
    dplyr::mutate(qc_status = dplyr::case_when(
      is.na(acquired_date) ~ "pending", qc_passed ~ "passed",
      .default = "flagged"
    )) %>%
    dplyr::select(-qc_passed)
  data.table::fwrite(run_log, file.path(folder, "run_log.csv"))
  folder
}

make_proteins <- function(n_proteins) {
  n_candidates <- 10
  genes <- c(
    established_neosubstrates,
    sprintf("SYNZF%02d", seq_len(n_candidates)),
    sprintf("SYNG%04d", seq_len(n_proteins - length(established_neosubstrates) - n_candidates))
  )
  is_neo <- genes %in% established_neosubstrates | grepl("^SYNZF", genes)
  tibble::tibble(
    uniprot_id = sprintf("S%05d", seq_along(genes)),
    gene_name = genes,
    scout_id = paste(gene_name, uniprot_id, sep = "_"),
    neosubstrate_category = dplyr::case_when(
      genes %in% established_neosubstrates ~ "established",
      grepl("^SYNZF", genes) ~ "candidate",
      .default = "none"
    ),
    zinc_finger = genes %in% zinc_finger_neosubstrates | grepl("^SYNZF", genes),
    # neosubstrates sit in the middle intensity range so their effects stay above the detection limit
    mean_log2_intensity = ifelse(is_neo, stats::rnorm(length(genes), 23, 0.7), stats::rnorm(length(genes), 22.5, 2)),
    noise_sd = stats::runif(length(genes), 0.08, 0.25)
  )
}

#' 96-well layout per plate: DMSO A1-F1, positive control G1-H1, 44 compounds in duplicate in columns 2-12 with the
#' replicates in row pairs (A/B, C/D, E/F, G/H). Measured column-wise.
make_plate_layout <- function() {
  purrr::map(seq_len(n_plates), \(p) {
    plate <- sprintf("P%d", p)
    controls <- tibble::tibble(
      compound = c(rep(vehicle_name, 6), rep(positive_control, 2)),
      replicate = c(paste0("R", 1:6), "R1", "R2"),
      row = 1:8, col = 1
    )
    k <- seq_len(compounds_per_plate)
    compounds <- tidyr::expand_grid(k = k, rep = 1:2) %>%
      dplyr::mutate(
        compound = sprintf("CMPD-%04d", (p - 1) * compounds_per_plate + k),
        replicate = paste0("R", rep),
        col = 2 + (k - 1) %/% 4,
        row = 2 * ((k - 1) %% 4) + rep
      ) %>%
      dplyr::select(compound, replicate, row, col)
    dplyr::bind_rows(controls, compounds) %>%
      dplyr::arrange(col, row) %>%
      dplyr::mutate(
        plate_id = plate,
        position = paste0(LETTERS[row], sprintf("%02d", col)),
        sample_id = sprintf("S%03d", (p - 1) * 96 + dplyr::row_number()),
        experiment_id = paste("MS01", "DEMO", "DIA", plate, sample_id, position, "DIA", "60SPD", "F1", "R1", sep = "_"),
        cell_line = "SYN-1 (synthetic)",
        compound_concentration = ifelse(compound == vehicle_name, "0.1 %", "1 uM"),
        time = "24 h"
      )
  }) %>%
    dplyr::bind_rows()
}

#' Run plan: one plate per week (plate 1 starts 13 days before today), plates beyond n_measured_plates still planned.
#' QC metrics and status are added once the plates are generated.
make_run_log <- function(layout, today) {
  start <- today - 13
  log <- layout %>%
    dplyr::mutate(plate_index = match(plate_id, unique(plate_id))) %>%
    dplyr::mutate(
      run_in_plate = dplyr::row_number(),
      planned_date = start + 7 * (plate_index - 1) + (run_in_plate - 1) %/% 20,
      .by = plate_id
    ) %>%
    dplyr::mutate(
      measured = plate_index <= n_measured_plates,
      acquired_date = dplyr::if_else(measured, planned_date + (plate_index == 2), as.Date(NA))
    )
  log %>%
    dplyr::transmute(plate_id, sample_id, position, compound, replicate, planned_date, acquired_date)
}

#' Simulated protein effects per compound of one plate (matrix proteins x compounds, log2 FC).
draw_effects <- function(proteins, compounds) {
  effect <- matrix(0, nrow(proteins), length(compounds), dimnames = list(proteins$gene_name, compounds))
  background <- which(proteins$neosubstrate_category == "none")
  for (cmpd in compounds) {
    profile <- compound_profiles[[cmpd]]
    if (is.null(profile)) {
      n_neo <- stats::rbinom(1, 4, 0.15)
      n_cand <- stats::rbinom(1, 1, 0.06)
      neo <- c(
        sample(names(neosubstrate_weights), n_neo, prob = neosubstrate_weights),
        sample(names(candidate_weights), n_cand, prob = candidate_weights)
      )
      profile <- list(
        neo = stats::setNames(-stats::runif(length(neo), 1.2, 3), neo),
        n_off = stats::rnbinom(1, size = 0.6, mu = 5)
      )
    }
    effect[names(profile$neo), cmpd] <- profile$neo
    off <- sample(background, profile$n_off)
    effect[off, cmpd] <- -stats::rgamma(length(off), shape = 4, rate = 3)
    up <- sample(setdiff(background, off), max(1, round(profile$n_off / 3)))
    effect[up, cmpd] <- stats::rgamma(length(up), shape = 4, rate = 3.5)
  }
  effect
}

#' Flat-header final table of one plate, in the column format prepare_final_table() parses, plus the QC call of
#' every run: passed when its protein IDs reach qc_min_share of the plate median. Flagged runs stay in the table.
#' @return list(final_table, qc = sample_id, qc_passed)
make_final_table <- function(proteins, plate_layout, underloaded, plate_offset, qc_min_share = 0.85) {
  samples <- plate_layout %>%
    dplyr::mutate(sample = paste(compound, replicate, sep = "|"))
  compounds <- setdiff(unique(samples$compound), vehicle_name)
  effect <- cbind(0, draw_effects(proteins, compounds))
  colnames(effect)[1] <- vehicle_name
  n_p <- nrow(proteins)
  n_s <- nrow(samples)

  sample_offset <- stats::rnorm(n_s, plate_offset, 0.12) - 2.5 * (samples$sample_id %in% underloaded)
  noise <- matrix(stats::rnorm(n_p * n_s), n_p, n_s) * proteins$noise_sd
  raw <- proteins$mean_log2_intensity + effect[, samples$compound] + noise + rep(sample_offset, each = n_p)
  raw[matrix(stats::runif(n_p * n_s), n_p, n_s) < stats::plogis(-(raw - 18.3) * 2) + 0.01] <- NA

  n_ids <- colSums(!is.na(raw))
  qc <- tibble::tibble(sample_id = samples$sample_id, qc_passed = n_ids >= qc_min_share * stats::median(n_ids))

  # median normalization, left-censored imputation of the treated samples (DMSO stays unimputed)
  normalized <- sweep(raw, 2, apply(raw, 2, stats::median, na.rm = TRUE) - 22)
  imputed <- normalized
  for (j in which(samples$compound != vehicle_name)) {
    miss <- is.na(normalized[, j])
    mu <- mean(normalized[, j], na.rm = TRUE) - 1.8 * stats::sd(normalized[, j], na.rm = TRUE)
    imputed[miss, j] <- stats::rnorm(sum(miss), mu, 0.3 * stats::sd(normalized[, j], na.rm = TRUE))
  }

  intensity_block <- function(type, m) {
    colnames(m) <- paste(type, samples$sample, sep = "|")
    tibble::as_tibble(round(2^m))
  }
  stats_prefix <- paste("Statistical testing results", "Compound vs DMSO", "Pairwise comparisons (t-test)", sep = "|")
  statistics <- purrr::map(compounds, \(cmpd) {
    moderated_t_test(
      imputed[, samples$compound == cmpd, drop = FALSE], normalized[, samples$compound == vehicle_name],
      is.na(normalized[, samples$compound == cmpd, drop = FALSE])
    ) %>%
      dplyr::rename_with(\(x) paste(stats_prefix, paste(cmpd, "-", vehicle_name), x, sep = "|"))
  }) %>%
    dplyr::bind_cols()

  final_table <- dplyr::bind_cols(
    tibble::tibble(
      `Identifier|Concatenated identifier` = proteins$scout_id,
      `Identifier|UniProt entry` = proteins$uniprot_id,
      `Identifier|Gene name` = proteins$gene_name,
      `Annotation|Plate ID` = plate_layout$plate_id[1],
      `Annotation|Neosubstrate category` = proteins$neosubstrate_category,
      `Annotation|Zinc finger` = proteins$zinc_finger
    ),
    intensity_block("Unnormalized LFQ intensities", raw),
    intensity_block("Normalized LFQ intensities", normalized),
    intensity_block("Imputed LFQ intensities", imputed),
    statistics
  )
  list(final_table = final_table, qc = qc)
}

#' Moderated t-test (pooled variance shrunk towards the median variance, prior df 4; in the spirit of limma's
#' empirical Bayes, Smyth 2004), BH-corrected. Stand-in for the pipeline's statistics, not a limma reimplementation.
#' Proteins with fewer than min_valid measured DMSO values are not tested.
#' @param treated matrix of imputed log2 intensities, @param vehicle matrix of normalized log2 intensities (with NA)
#' @param treated_imputed logical matrix, TRUE where a treated value was imputed
moderated_t_test <- function(treated, vehicle, treated_imputed, min_valid = 3, prior_df = 4) {
  n_t <- ncol(treated)
  n_v <- rowSums(!is.na(vehicle))
  ss <- rowSums((treated - rowMeans(treated))^2) + rowSums((vehicle - rowMeans(vehicle, na.rm = TRUE))^2, na.rm = TRUE)
  df <- n_t + n_v - 2
  s2 <- ss / df
  s2_prior <- stats::median(s2[n_v >= min_valid], na.rm = TRUE)
  s2_post <- (prior_df * s2_prior + df * s2) / (prior_df + df)
  log2fc <- rowMeans(treated) - rowMeans(vehicle, na.rm = TRUE)
  t <- log2fc / sqrt(s2_post * (1 / n_t + 1 / n_v))
  p <- 2 * stats::pt(-abs(t), prior_df + df)
  tested <- n_v >= min_valid
  p[!tested] <- NA
  valid <- paste(
    apply(ifelse(treated_imputed, "i", "v"), 1, paste, collapse = ""),
    apply(ifelse(is.na(vehicle), "m", "v"), 1, paste, collapse = ""),
    sep = " - "
  )
  tibble::tibble(
    `log2 fold change` = signif(ifelse(tested, log2fc, NA), 5),
    `t-statistic` = signif(ifelse(tested, t, NA), 5),
    `p-value (uncorrected)` = signif(p, 5),
    `p-value (corrected)` = signif(stats::p.adjust(p, method = "BH"), 5),
    `Valid value string` = valid
  )
}
