# infer_plate_layout() copied unchanged from the qmd-reports template (R/fcn_plate_view.R),
# without the ggplot2 plate plots of that file. Keep in sync with the template.
require(magrittr)


#' Infer the plate layout from the raw file in the sample annotation
#'
#' @param sample_annotation data frame with sample annotation
#' @param project_params list with project parameters
#' @param measurement_order_col name column that indicates the measurement order, default to NULL assuming that the order was column wise on the plate
#' @param columns_to_keep list of columns that should be kept from the sample annotation (outside the group headers), defaults to c()
#' @returns data frame similar to the sample annotation but with new columns that indicate the position of the sample on the plate
#'
#' @examples
#' plate <- infer_plate_layout(sample_annotation, project_params, measurement_order_col = "measurement_order",
#'                             columns_to_keep = "cell_line")
infer_plate_layout <- function(sample_annotation, project_params, measurement_order_col = NULL, 
                               columns_to_keep = c()) {
  
  tryCatch(expr = {
    plate_assignment <- sample_annotation %>%
      janitor::clean_names() %>%
      tidyr::separate(experiment_id, into = c("machine", "project_id", "workflow", "plate_id", "sample_id", "position", 
                                              "ms_method_1", "ms_method_2", "fraction", "experiment_replicate"), sep = "_") %>%
      tidyr::unite(ms_method, c("ms_method_1", "ms_method_2")) %>%
      tidyr::unite(group, c(project_params$group_header, project_params$additional_header), remove = FALSE) %>%
      tidyr::unite(sample, c("group", "replicate"), remove = FALSE) %>%
      dplyr::mutate_at(project_params$additional_header, \(x) as.numeric(gsub(project_params$additional_header_unit, "", x))) %>%
      dplyr::select(any_of(c("sample", "group", columns_to_keep, project_params$group_header, project_params$additional_header, 
                             "replicate", "position", measurement_order_col)))
    
    well_positions <- expand.grid(letter = LETTERS[c(1:8)], idx = c(1:12)) %>%
      dplyr::mutate(letter_idx = purrr::map_int(letter, \(x) grep(x, LETTERS)),
                    position = purrr::map2_chr(letter, idx, \(x, n) paste0(x, stringr::str_pad(n, 2, pad = "0"))))
    
    plate_assignment <- plate_assignment %>%
      dplyr::left_join(well_positions, by = "position")
    
    if (is.null(measurement_order_col) | isTRUE(!measurement_order_col %in% names(plate_assignment))) {
      message("Measurement order was inferred from plate layout (column-wise).")
      plate_assignment <- plate_assignment %>%
        dplyr::arrange(idx, letter_idx) %>%
        dplyr::mutate(measurement_order = 1:nrow(plate_assignment))
      measurement_order_tag <- "column-wise"
    } else {
      message(paste0("Measurement order was adopted from sample annotation (", measurement_order_col, ")."))
      plate_assignment <- plate_assignment %>%
        dplyr::rename(measurement_order = .data[[measurement_order_col]])
      measurement_order_tag <- "from sample annotation"
    }
  }, 
  error = function (e) {
    message(paste(e, "Sample annotation could not be parsed.\n", sep = "\n"))
    plate_assignment <- NULL
    measurement_order_tag <- NULL
  })
  
  return(list(plate_assignment = plate_assignment, measurement_order_tag = measurement_order_tag))
}
