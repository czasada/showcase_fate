# Numbers behind the Experimental setup tab: one row per well of every plate (planned plates included) and the
# study design, taken from the sample annotation (via plate_assignment) and the run log.

#' All wells of all plates: run log joined with the measured plate layout.
setup_wells <- function(study) {
  layout <- study$plate_assignment %>%
    dplyr::select(plate_id, position, sample, measurement_order)
  study$run_log %>%
    dplyr::left_join(layout, by = c("plate_id", "position")) %>%
    dplyr::mutate(
      sample = dplyr::coalesce(sample, paste(compound, replicate, sep = "_")),
      plate_row = substr(position, 1, 1),
      plate_col = as.integer(substr(position, 2, 3)),
      well_type = dplyr::case_when(
        compound == vehicle_name ~ "DMSO",
        compound == positive_control ~ "positive control",
        .default = "compound"
      ),
      run_status = dplyr::case_when(
        qc_status == "flagged" ~ "flagged (below QC threshold)",
        qc_status == "pending" ~ "planned",
        .default = qc_status
      )
    )
}

#' Study design rows (label, value) from the sample annotation and the run log.
design_summary <- function(study) {
  wells <- setup_wells(study)
  annotation <- study$plate_assignment
  first <- function(x) paste(unique(x), collapse = ", ")
  per_plate <- wells %>%
    dplyr::count(plate_id, well_type) %>%
    dplyr::summarise(n = paste(unique(n), collapse = "/"), .by = well_type)
  per_plate_n <- function(type) per_plate$n[per_plate$well_type == type]
  n_replicates <- wells %>%
    dplyr::filter(well_type == "compound") %>%
    dplyr::count(plate_id, compound) %>%
    dplyr::pull(n) %>%
    unique()
  treated <- annotation$compound != vehicle_name
  tibble::tribble(
    ~label, ~value,
    "Cell line", first(annotation$cell_line),
    "Treatment", sprintf(
      "%s, %s (DMSO %s)", first(annotation$compound_concentration[treated]), first(annotation$time),
      first(annotation$compound_concentration[!treated])
    ),
    "Plates", sprintf("%d × 96 wells", dplyr::n_distinct(wells$plate_id)),
    "Per plate", sprintf(
      "%s DMSO · %s positive control · %d compounds × %s", per_plate_n("DMSO"), per_plate_n("positive control"),
      as.integer(per_plate_n("compound")) / max(n_replicates), paste(n_replicates, collapse = "/")
    ),
    "Compounds", as.character(dplyr::n_distinct(wells$compound[wells$well_type == "compound"])),
    "Reference", "DMSO wells of the same plate",
    "Measurement", sprintf("%s, label-free, column-wise", first(sub("_", " ", annotation$ms_method)))
  )
}

#' Frames around the two wells of each compound (plotly shapes on the plate map's category y axis).
replicate_pair_shapes <- function(wells, rows) {
  pairs <- wells %>%
    dplyr::filter(well_type == "compound") %>%
    dplyr::mutate(y = match(plate_row, rows) - 1) %>%
    dplyr::summarise(y0 = min(y) - 0.45, y1 = max(y) + 0.45, x = dplyr::first(plate_col), .by = compound)
  purrr::pmap(pairs, \(compound, y0, y1, x) {
    list(
      type = "rect", x0 = x - 0.45, x1 = x + 0.45, y0 = y0, y1 = y1, yref = "y",
      line = list(color = "#6c757d", width = 1), fillcolor = "rgba(0,0,0,0)"
    )
  })
}
