# Synthetic degradomics data for the public demo dashboard.
# Everything here is generated from a fixed seed - no project data is read. Structure mimics the
# ProteomeScout final table (intensities, statistics, annotation) and the WP1 degradomics report.

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
deep_dive_compounds <- c("CMPD-0042", "CMPD-0117", "CMPD-0158", "CMPD-0203", "CMPD-0266")

# Compound profiles: named log2 FC per neosubstrate, plus how many random off-targets go down/up.
compound_profiles <- list(
  "CMPD-0042" = list(neo = c(IKZF1 = -2.6, IKZF3 = -2.9, ZFP91 = -1.4, ZNF692 = -1.1), n_off_down = 18, n_off_up = 6),
  "CMPD-0117" = list(neo = c(GSPT1 = -2.3, GSPT2 = -1.5), n_off_down = 35, n_off_up = 20),
  "CMPD-0158" = list(neo = c(CSNK1A1 = -2.1, IKZF1 = -1.6, IKZF3 = -1.2), n_off_down = 12, n_off_up = 5),
  "CMPD-0203" = list(
    neo = c(
      IKZF1 = -2.2, IKZF2 = -1.9, IKZF3 = -2.4, ZFP91 = -1.7, ZNF692 = -1.5, SALL4 = -1.3,
      WIZ = -1.8, FIZ1 = -1.2, ZNF827 = -1.1, SYNZF03 = -1.6, SYNZF07 = -1.3
    ),
    n_off_down = 25, n_off_up = 8
  ),
  "CMPD-0266" = list(neo = c(PDE6D = -0.8), n_off_down = 4, n_off_up = 3)
)

#' Build the full synthetic data set: protein annotation, sample annotation, intensities,
#' statistics (each compound vs DMSO) and the 384-well screen status table.
make_demo_data <- function(seed = 640, n_proteins = 5000, n_replicates = 4) {
  set.seed(seed)

  proteins <- make_protein_annotation(n_proteins)
  samples <- make_sample_annotation(n_replicates)
  intensities <- make_intensities(proteins, samples)
  statistics <- compute_statistics(intensities, samples)
  screen <- make_screen(statistics, intensities, proteins)

  list(proteins = proteins, samples = samples, intensities = intensities, statistics = statistics, screen = screen)
}

make_protein_annotation <- function(n_proteins) {
  n_candidates <- 10
  genes <- c(
    established_neosubstrates,
    sprintf("SYNZF%02d", seq_len(n_candidates)),
    sprintf("SYNG%04d", seq_len(n_proteins - length(established_neosubstrates) - n_candidates))
  )
  tibble::tibble(
    protein_id = sprintf("SYN%05d", seq_along(genes)),
    gene_name = genes,
    neosubstrate_category = dplyr::case_when(
      genes %in% established_neosubstrates ~ "established",
      grepl("^SYNZF", genes) ~ "candidate",
      .default = "none"
    ),
    zinc_finger = genes %in% zinc_finger_neosubstrates | grepl("^SYNZF", genes),
    # protein-level abundance and noise; neosubstrates sit in the lower-middle intensity range
    mean_log2_intensity = ifelse(genes %in% c(established_neosubstrates, sprintf("SYNZF%02d", seq_len(n_candidates))),
      stats::rnorm(length(genes), 23, 0.7), stats::rnorm(length(genes), 22.5, 2)
    ),
    noise_sd = stats::runif(length(genes), 0.08, 0.25)
  )
}

make_sample_annotation <- function(n_replicates) {
  conditions <- c(vehicle_name, deep_dive_compounds)
  wells <- expand.grid(row = LETTERS[1:8], col = 1:12, stringsAsFactors = FALSE)
  samples <- tidyr::expand_grid(compound = conditions, replicate = seq_len(n_replicates)) %>%
    dplyr::mutate(
      sample = paste0(compound, "_R", replicate),
      cell_line = "SYN-1 (synthetic)",
      compound_concentration = ifelse(compound == vehicle_name, "0.1 %", "1 uM"),
      time = "24 h",
      # replicates in rows, conditions in columns 2..7, like a 96-well plate run
      plate_row = LETTERS[replicate],
      plate_col = match(compound, conditions) + 1,
      well = paste0(plate_row, plate_col),
      exclude = FALSE
    )
  samples
}

make_intensities <- function(proteins, samples) {
  n_p <- nrow(proteins)
  n_s <- nrow(samples)

  # true log2 fold changes per protein and condition
  effect <- matrix(0, n_p, n_s, dimnames = list(proteins$protein_id, samples$sample))
  for (cmpd in names(compound_profiles)) {
    profile <- compound_profiles[[cmpd]]
    cols <- samples$compound == cmpd
    neo_idx <- match(names(profile$neo), proteins$gene_name)
    effect[neo_idx, cols] <- profile$neo
    background <- which(proteins$neosubstrate_category == "none")
    off_down <- sample(background, profile$n_off_down)
    off_up <- sample(setdiff(background, off_down), profile$n_off_up)
    effect[off_down, cols] <- -stats::rgamma(length(off_down), shape = 4, rate = 3.5)
    effect[off_up, cols] <- stats::rgamma(length(off_up), shape = 4, rate = 4)
  }

  sample_offset <- stats::rnorm(n_s, 0, 0.12)
  sample_offset[samples$sample == paste0(deep_dive_compounds[4], "_R3")] <- -0.7 # one under-loaded sample
  noise <- matrix(stats::rnorm(n_p * n_s), n_p, n_s) * proteins$noise_sd
  raw <- proteins$mean_log2_intensity + effect + noise + rep(sample_offset, each = n_p)

  # left-censored missingness (low intensities drop out) plus a little random missingness
  p_missing <- stats::plogis(-(raw - 18.3) * 2) + 0.01
  raw[matrix(stats::runif(n_p * n_s), n_p, n_s) < p_missing] <- NA

  # median normalization to the global median
  sample_medians <- apply(raw, 2, stats::median, na.rm = TRUE)
  normalized <- sweep(raw, 2, sample_medians - mean(sample_medians))

  # left-censored imputation per sample (down-shifted normal distribution)
  imputed <- normalized
  for (j in seq_len(n_s)) {
    miss <- is.na(normalized[, j])
    mu <- mean(normalized[, j], na.rm = TRUE) - 1.8 * stats::sd(normalized[, j], na.rm = TRUE)
    imputed[miss, j] <- stats::rnorm(sum(miss), mu, 0.3 * stats::sd(normalized[, j], na.rm = TRUE))
  }

  tibble::tibble(
    protein_id = rep(proteins$protein_id, times = n_s),
    sample = rep(samples$sample, each = n_p),
    raw_log2_intensity = as.vector(raw),
    normalized_log2_intensity = as.vector(normalized),
    imputed_log2_intensity = as.vector(imputed),
    is_imputed = is.na(as.vector(normalized))
  ) %>%
    dplyr::left_join(samples %>% dplyr::select(sample, compound, replicate), by = "sample")
}

#' Welch t-test of each compound vs DMSO on imputed log2 intensities, BH-corrected per comparison.
#' A protein is tested when it has at least 3 measured (not imputed) values in either group.
compute_statistics <- function(intensities, samples, min_valid = 3) {
  wide <- function(col) {
    m <- tidyr::pivot_wider(intensities, id_cols = protein_id, names_from = sample, values_from = dplyr::all_of(col))
    as.matrix(tibble::column_to_rownames(m, "protein_id"))
  }
  values <- wide("imputed_log2_intensity")
  measured <- !wide("is_imputed")
  vehicle_cols <- samples$sample[samples$compound == vehicle_name]
  compounds <- setdiff(unique(samples$compound), vehicle_name)

  purrr::map(compounds, \(cmpd) {
    treated_cols <- samples$sample[samples$compound == cmpd]
    v <- values[, vehicle_cols, drop = FALSE]
    t <- values[, treated_cols, drop = FALSE]
    n_v <- ncol(v)
    n_t <- ncol(t)
    var_v <- row_var(v)
    var_t <- row_var(t)
    valid_v <- rowSums(measured[, vehicle_cols, drop = FALSE])
    valid_t <- rowSums(measured[, treated_cols, drop = FALSE])
    tibble::tibble(
      protein_id = rownames(values),
      comparison = paste(cmpd, "-", vehicle_name),
      compound = cmpd,
      log2_fold_change = rowMeans(t) - rowMeans(v),
      se = sqrt(var_t / n_t + var_v / n_v),
      df = se^4 / ((var_t / n_t)^2 / (n_t - 1) + (var_v / n_v)^2 / (n_v - 1)),
      n_imputed = (n_t - valid_t) + (n_v - valid_v),
      tested = valid_t >= min_valid | valid_v >= min_valid
    ) %>%
      dplyr::filter(tested) %>%
      dplyr::mutate(
        p_value_uncorrected = 2 * stats::pt(-abs(log2_fold_change / se), df),
        p_value_corrected = stats::p.adjust(p_value_uncorrected, method = "BH")
      ) %>%
      dplyr::select(
        protein_id, comparison, compound, log2_fold_change, p_value_uncorrected, p_value_corrected, n_imputed
      )
  }) %>%
    dplyr::bind_rows()
}

row_var <- function(m) {
  rowSums((m - rowMeans(m))^2) / (ncol(m) - 1)
}

screen_stages <- c("planned", "acquired", "QC failed", "QC passed", "analysed", "reported")

#' 384-well screen: DMSO in column 1, positive control / blanks in column 24, 352 compounds in 2-23.
#' Acquisition proceeds column by column, so earlier columns are further along. Deep-dive compounds sit
#' in reported wells and take their numbers from the detailed statistics; all other numbers are drawn.
make_screen <- function(statistics, intensities, proteins, screen_hit_p = 0.05, screen_hit_log2fc = -1) {
  wells <- expand.grid(plate_row = LETTERS[1:16], plate_col = 1:24, stringsAsFactors = FALSE) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(plate_col, plate_row) %>%
    dplyr::mutate(well = paste0(plate_row, plate_col))

  compound_wells <- wells$plate_col %in% 2:23
  wells$well_type <- dplyr::case_when(
    wells$plate_col == 1 ~ "vehicle",
    wells$plate_col == 24 & wells$plate_row %in% LETTERS[1:8] ~ "positive control",
    wells$plate_col == 24 ~ "blank",
    .default = "compound"
  )
  wells$compound <- dplyr::case_when(
    wells$well_type == "vehicle" ~ vehicle_name,
    wells$well_type == "positive control" ~ "POS-CTRL",
    wells$well_type == "blank" ~ "blank",
    .default = NA_character_
  )
  ids <- sprintf("CMPD-%04d", seq_len(sum(compound_wells)))
  wells$compound[compound_wells] <- ids

  # deep-dive compounds go to fixed wells in reported columns (swap with the compound sitting there)
  deep_wells <- c("C3", "H4", "K6", "E7", "N9")
  for (i in seq_along(deep_dive_compounds)) {
    swap <- c(which(wells$compound == deep_dive_compounds[i]), which(wells$well == deep_wells[i]))
    wells$compound[swap] <- wells$compound[rev(swap)]
  }

  # stage by column: 2-9 reported, 10-13 analysed, 14-17 QC passed, 18-20 acquired, 21-23 planned
  col_stage <- cut(
    wells$plate_col, c(0, 9, 13, 17, 20, 24),
    labels = c("reported", "analysed", "QC passed", "acquired", "planned")
  )
  wells$stage <- as.character(col_stage)
  wells$stage[wells$well_type %in% c("vehicle", "positive control", "blank")] <- "reported"
  acquired_plus <- wells$stage != "planned"
  qc_fail <- acquired_plus & wells$well_type == "compound" & stats::runif(nrow(wells)) < 0.03 &
    !wells$compound %in% deep_dive_compounds
  wells$stage[qc_fail] <- "QC failed"
  wells$stage[wells$well_type == "blank"] <- "acquired"

  n <- nrow(wells)
  wells$n_ids <- ifelse(acquired_plus, round(stats::rnorm(n, 4650, 90)), NA)
  wells$n_ids[qc_fail] <- round(stats::rnorm(sum(qc_fail), 3100, 300))
  wells$n_ids[wells$well_type == "blank"] <- round(stats::runif(sum(wells$well_type == "blank"), 20, 120))

  analysed <- wells$stage %in% c("analysed", "reported") & wells$well_type %in% c("compound", "positive control")
  wells$n_degraded <- NA_integer_
  wells$n_degraded[analysed] <- stats::rnbinom(sum(analysed), size = 0.6, mu = 6)
  wells$n_degraded[wells$well_type == "positive control"] <- stats::rpois(8, 20)
  wells$n_neosubstrates_degraded <- NA_integer_
  n_neo <- stats::rbinom(sum(analysed), 4, 0.12) * (wells$n_degraded[analysed] > 0)
  wells$n_neosubstrates_degraded[analysed] <- pmin(n_neo, wells$n_degraded[analysed])
  wells$n_neosubstrates_degraded[wells$well_type == "positive control"] <- 3
  wells$top_neosubstrate <- NA_character_
  has_neo <- analysed & wells$n_neosubstrates_degraded > 0
  wells$top_neosubstrate[has_neo] <- sample(established_neosubstrates[1:12], sum(has_neo), replace = TRUE)
  wells$top_neosubstrate[wells$well_type == "positive control"] <- "IKZF1"

  # deep-dive compounds: use their actual numbers
  hits <- statistics %>%
    dplyr::filter(p_value_corrected <= screen_hit_p, log2_fold_change <= screen_hit_log2fc) %>%
    dplyr::left_join(proteins %>% dplyr::select(protein_id, gene_name, neosubstrate_category), by = "protein_id")
  for (i in seq_along(deep_dive_compounds)) {
    cmpd_hits <- hits %>% dplyr::filter(compound == deep_dive_compounds[i])
    n_ids_cmpd <- intensities %>%
      dplyr::filter(compound == deep_dive_compounds[i], !is_imputed) %>%
      dplyr::distinct(protein_id) %>%
      nrow()
    row <- which(wells$well == deep_wells[i])
    wells$n_ids[row] <- n_ids_cmpd
    wells$n_degraded[row] <- nrow(cmpd_hits)
    wells$n_neosubstrates_degraded[row] <- sum(cmpd_hits$neosubstrate_category == "established")
    wells$top_neosubstrate[row] <- cmpd_hits %>%
      dplyr::filter(neosubstrate_category == "established") %>%
      dplyr::slice_min(log2_fold_change, n = 1) %>%
      dplyr::pull(gene_name) %>%
      dplyr::first()
  }

  wells %>%
    dplyr::mutate(
      stage = factor(stage, levels = screen_stages),
      deep_dive = compound %in% deep_dive_compounds
    )
}
