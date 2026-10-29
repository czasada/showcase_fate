# Read the study the way the reports do: one final table per plate via prepare_final_table(), the plate layout via
# infer_plate_layout(), plus the run log. The dashboard only aggregates these tables; it runs no analysis itself.

# Project parameters as parse_project_config() would return them for a compound x replicate design.
demo_project_params <- function(final_table_path) {
  list(
    header = c("compound", "replicate"), group_header = "compound", additional_header = NULL,
    additional_header_unit = NULL, quant_type = "dia_lfq", final_table_path = final_table_path
  )
}

#' Read all plates in a study folder and stack them.
#' @return list(identifier, annotations, intensities, statistics, plate_assignment, run_log); intensities, statistics
#'   and plate_assignment carry plate_id, statistics also compound (group 1 of the comparison)
read_study <- function(folder) {
  paths <- list.files(folder, pattern = "^final_table_.*_flat_header\\.csv$", full.names = TRUE)
  plates <- purrr::map(paths, \(path) {
    params <- demo_project_params(path)
    data <- suppressMessages(prepare_final_table(params))
    plate <- unique(data$annotations$plate_id)
    stopifnot(length(plate) == 1)
    annotation_path <- file.path(folder, paste0("sample_annotation_", plate, ".csv"))
    layout <- suppressMessages(infer_plate_layout(
      data.table::fread(annotation_path), params,
      columns_to_keep = c("plate_id", "cell_line", "compound_concentration", "time", "ms_method")
    ))
    list(
      identifier = data$identifier,
      annotations = data$annotations,
      intensities = dplyr::mutate(data$intensities, plate_id = plate),
      statistics = data$statistics_t %>%
        dplyr::mutate(plate_id = plate, compound = sub(" - .*$", "", comparison)),
      plate_assignment = layout$plate_assignment
    )
  })
  stack <- function(name) dplyr::bind_rows(purrr::map(plates, name))
  list(
    identifier = dplyr::distinct(stack("identifier")),
    annotations = stack("annotations") %>%
      dplyr::select(-plate_id) %>%
      dplyr::distinct() %>%
      dplyr::mutate(zinc_finger = as.logical(zinc_finger)),
    intensities = stack("intensities"),
    statistics = stack("statistics"),
    plate_assignment = stack("plate_assignment"),
    run_log = data.table::fread(file.path(folder, "run_log.csv")) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(dplyr::across(dplyr::ends_with("_date"), as.Date))
  )
}

#' Synthetic study written to a temporary folder and read back.
load_demo_study <- function(folder = file.path(tempdir(), "demo_study"), ...) {
  write_synthetic_study(folder, ...)
  read_study(folder)
}
