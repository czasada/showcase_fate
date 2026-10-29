## FUNCTIONS TO READ AND PREPARE THE PROJECT FINAL TABLE
# Copied unchanged from the qmd-reports template. Keep in sync.
require(magrittr)


#' Add the final table path to the project parameter list from the parameters of the qmd file, i.e. final_table_tag and
#' final_table_path
#' 
#' If the latter is NULL the path is inferred from the project parameters
#'
#' @param project_params list of project parameters, output of parse_project_config()
#' @param qmd_params list of input parameters of the qmd file, can contain a specific path to the final table, 
#' if not the path is deduced from the project parameters (analysis_folder + 06_final_tables)
#' @return updated project paramter list with the path to the final table as final_table_path
#' 
#' @examples
#' project_params <- parse_project_config(project_config) %>%
#'   add_final_table_path(params)
add_final_table_path <- function(project_params, qmd_params) {
  
  if (is.null(qmd_params$final_table_path)) {
    final_table_folder <- file.path(project_params$analysis_folder, "06_final_tables")
    final_table_file <- list.files(final_table_folder, 
                                   pattern = paste0("(.*)", qmd_params$final_table_tag, "(.*)_flat_header.csv"))
    final_table_path <- file.path(final_table_folder, final_table_file)
  } else {
    final_table_path <- file.path(qmd_params$final_table_path)
  }
  
  if (length(final_table_path) > 1) {
    warning(paste0("\nThe final table is not uniquely identified. Possible tables are:\n\t", paste(final_table_path, collapse = "\n\t")))
  }
  message(paste0("\nThe final table path that is added to the project parameters:\n\t", final_table_path[1], ".\n",
                 "If this is wrong specify final_table_tag or the correct full final_table_path in the input parameters.\n"))
  
  project_params$final_table_tag <- qmd_params$final_table_tag
  project_params$final_table_path <- final_table_path[1]
  
  return(project_params)
}


#' Read in the final table given its file path
#'
#' @param final_table_path path to flat header final table
#' @return final table as in csv file
#' 
#' @examples
#' final_table <- read_final_table(project_params$final_table_path)
read_final_table <- function(final_table_path) {
  
  final_table <- data.table::fread(final_table_path)
  colnames(final_table) <- gsub("[\\|]+", "|", colnames(final_table))
  
  return(final_table)
}


#' Rename the ID columns of the final table
#'
#' @param final_table final table 
#' @return updated final table with renamed identifier columns
#'
#' @examples
#' final_table <- read_final_table(project_params, params) %>%
#'   prepare_id_columns(project_params)
prepare_id_columns <- function (final_table) {
  
  final_table_columns_renaming <- c(scout_id = "Identifier|Concatenated identifier",
                                    uniprot_id = "Identifier|UniProt entry",
                                    uniprot_entry_name = "Identifier|UniProt entry name",
                                    gene_name = "Identifier|Gene name",
                                    mod_aa = "Identifier|Modification associated amino acid(s)",
                                    mod_position = "Identifier|Modification associated amino acid(s)",
                                    mod_position = "Identifier|Modification site position(s) (UniProt coordinates)",
                                    pep_sequence_with_mod = "Identifier|Peptide sequence with modification")
  
  final_table_with_id_cols <- final_table %>%
    dplyr::rename(any_of(final_table_columns_renaming)) %>%
    dplyr::relocate(scout_id) %>%
    dplyr::filter(scout_id != "") %>%
    tidyr::separate(uniprot_id, c("uniprot_id_noiso", "iso_n"), sep = "-", remove = FALSE, 
                    extra = "drop", fill = "right") %>%
    tidyr::separate(scout_id, c("protein_scout_id", "pep_sequence"), sep = "\\|", remove = FALSE,
                    extra = "drop", fill = "right") %>%
    janitor::remove_empty(which = "cols")
  
  if (nrow(final_table_with_id_cols) < nrow(final_table)) {
    n_diff <- nrow(final_table) - nrow(final_table_with_id_cols)
    warning(paste("There are", n_diff, "rows in the final table with empty Concatenated identifier that are filtered out!"))
  }
  
  return(final_table_with_id_cols)
}


#' Extract the identifier columns from a table. Identifier columns can be "scout_id", "uniprot_id", "uniprot_id_noiso", 
#' "iso_n", "uniprot_entry_name", "gene_name", "mod_aa", "mod_position", "pep_sequence"or "pep_sequence_with_mod"
#'
#' @param table any table
#' @return intersection between the column names of the table and the theoretical identifier named above
#'
#' @examples
#' id_columns <- get_id_columns(final_table %>% prepare_id_columns(project_params))
get_id_columns <- function(table) {
  
  theoretical_ids <- c("scout_id", "uniprot_id", "uniprot_id_noiso", "iso_n", "uniprot_entry_name", 
                       "gene_name", "mod_aa", "mod_position", "pep_sequence", "pep_sequence_with_mod")            
  id_columns <- dplyr::intersect(theoretical_ids, colnames(table))
  
  return(id_columns)
}


#' Prepare part of the final table as a smaller long-format table with the levels of the final table (separated by "|")
#' split into separate columns
#' 
#' If wanted the last level is kept wide, and certain column levels are dropped
#'
#' @param final_table final table 
#' @param columns_to_pick column names that should be part of the subtable besides scout_id
#' @param column_levels list with the description of the levels into which the column header are separated
#' (separator = "|")
#' @param leave_last_level_wide TRUE if the last part of the column template should stay in wide format, FALSE if the
#' table should be completely long with the values in a column called "value", defaults to FALSE
#' @param drop_levels list with elements of column_levels that should be dropped while processing, defaults to c()
#' @return subtable of final table in long format with separated column levels
#'
#' @examples
#' intensities <- prepare_final_table_section(final_table, 
#'                                            columns_to_pick = c("Normalized TMT reporter intensities|plex01_PP|Cmpd1|1hr|R1",
#'                                                                "Normalized TMT reporter intensities|plex01_PP|Cmpd1|1hr|R2",
#'                                                                "Normalized TMT reporter intensities|plex01_PP|Cmpd1|1hr|R3"),
#'                                            column_levels = c("type", "plex", "compound", "time", "replicate"),
#'                                            leave_last_level_wide = FALSE, drop_levels = c())
prepare_final_table_section <- function (final_table, columns_to_pick, column_levels, leave_last_level_wide = FALSE,
                                         drop_levels = c()) {
  
  subtable <- final_table %>%
    dplyr::select(scout_id, all_of(columns_to_pick))
  
  idx_to_keep <- which(!column_levels %in% drop_levels)
  column_levels <- setdiff(column_levels, drop_levels)
  colnames(subtable) <- purrr::map_chr(colnames(subtable), \(x) {
    if (x == "scout_id") {x} else {paste(stringr::str_split_1(x, pattern = "\\|")[idx_to_keep], collapse = "|")}
  })
  
  if (length(unique(colnames(subtable))) != length(colnames(subtable))) {
    stop(paste0("The subtable contains duplicate column names after dropping c(", paste(drop_levels, collapse = ", "), ")"))
  }
  
  if (length(column_levels) <= 1) {
    return (subtable %>% janitor::clean_names())
  } else {
    if (leave_last_level_wide == TRUE) {
      subtable_long <- subtable %>% 
        tidyr::pivot_longer(-scout_id, names_to = c("column_template", ".value"), names_sep = "\\|\\s*(?=[^\\|]+$)") %>%
        dplyr::filter_at(dplyr::vars(-scout_id, -column_template), dplyr::any_vars(!is.na(.) & !is.infinite(.)))
      column_levels <- head(column_levels, n = length(column_levels) - 1)
    } else {
      subtable_long <- subtable %>%
        tidyr::pivot_longer(-scout_id, names_to = "column_template", values_to = "value") %>%
        dplyr::filter(!is.na(value) & !is.infinite(value))
    }
    subtable_long_split <- subtable_long %>%
      tidyr::separate_wider_delim(column_template, delim = "|", names = column_levels)
    return(subtable_long_split %>% janitor::clean_names())
  }
}


#' Reads in the final table and creates a named list, dividing the final table into subtables by
#' (1) identifier
#' (2) annotations
#' (3) msms and peptide counts, peptide sequences and coverage (if present in final table)
#' (4) intensities (contains unnormalized, normalized and imputed intensities if present in final table)
#' (5) scaled intensities 
#' (6) curve parameters
#' (7) CR parameters
#' (8) statistics (t-test)
#' (9) statistics (F-test)
#' (10) F-test for curve comparison
#' 
#' If a subtable does not exist the respective element in the list is NULL
#' 
#' Sample groups are defined using the header levels of the project parameters (and heavily depend on the group header
#' and their order)
#'
#' @param project_params list of project parameters, output of parse_project_config() with add_final_table_path(), 
#' important elements are: final_table_path, group_header, additional_header and quant_type
#' @param keep_group_header if TRUE the group header columns are joined to a new column 'group' but are not removed 
#' from the tables (affects intensities, cr and curve parameters), if FALSE only the combined column 'group' remains in 
#' the output tables, defaults to TRUE
#' @return a list of subtables, all 10 elements are present, but if there is no info in the final table or it couldn't 
#' be processed, the subtable is an empty tibble
#'
#' @examples
#' data <- prepare_final_table(project_params$final_table_path)
prepare_final_table <- function(project_params, keep_group_header = TRUE) {
  
  final_table <- read_final_table(project_params$final_table_path) %>%
    prepare_id_columns()
  final_table_columns <- colnames(final_table)
  extra_level_info <- if (grepl("tmt", project_params$quant_type)) "plex" else "int_type"
  replicate_info <- if ("replicate" %in% project_params$header) "replicate" else NULL
  
  # (1) IDENTIFIER
  id_columns <- get_id_columns(final_table)
  identifier <- final_table %>%
    dplyr::select(all_of(id_columns))
  
  # (2) ANNOTATIONS
  # unfortunately need some clean up..
  # the manual annotations can be under Identifier but also under Annotation..
  annotations_regex <- paste0("^Annotation(.*)|(.*)Protein of interest(.*)|",
                              "^Identification\\/quantification reliability(.*)Class 1 site(.*)|",
                              "^Identifier(.*)([Aa]nnotation|[Pp]athway|[Pp]rotease|[Kk]inase|[Kk]inome|ATP binding|Direct binder|Might be coenriched)(.*)")
  annotations_columns <- regmatches(final_table_columns, regexpr(annotations_regex, final_table_columns))
  annotations <- final_table %>%
    dplyr::select(scout_id, all_of(annotations_columns)) %>%
    dplyr::rename_with(~ stringr::str_replace(., "Identifier\\|", "")) %>%
    dplyr::rename_with(~ stringr::str_replace(., "Identification\\/quantification reliability\\|", "")) %>%
    dplyr::rename_with(~ stringr::str_replace(., "Annotation\\|", "")) %>%
    dplyr::rename_with(~ stringr::str_replace(., "PhosphoSitePlus\\|PSP", "PSP")) %>%
    dplyr::rename_with(~ stringr::str_replace(., "PhosphoSitePlus\\|", "PSP")) %>%
    dplyr::rename_with(~ stringr::str_replace(., "Gene Ontology\\|Associated with ", "")) %>%
    dplyr::rename_with(~ stringr::str_replace(., "Pathways\\|Member of ", "")) %>%
    dplyr::rename_with(~ stringr::str_replace(., "Proteins of interest\\|", "")) %>%
    dplyr::rename_with(~ stringr::str_replace(., "Additional annotations\\|", "")) %>%
    janitor::clean_names()
  
  # (3) PEPTIDE AND MSMS COUNT
  # is rather inconsistent across scout types/experimental designs..
  peptide_metrices <- "MS\\s?[Cc]ount|[Pp]eptides|Peptide sequences|peptide\\_sequences|[Ss]equence coverage \\[\\%\\]|Missed cleavages"
  peptide_count_regex <- paste0("^Identification/quantification reliability(.*)(", peptide_metrices, ").*")
  peptide_count_columns <- regmatches(final_table_columns, regexpr(peptide_count_regex, final_table_columns))

  peptide_counts <- tryCatch(
    expr = { 
      n_levels <- unique(purrr::map_int(peptide_count_columns, \(x) length(stringr::str_split_1(x, "\\|"))))
      peptide_level_last <- all(stringr::str_ends(peptide_count_columns, pattern = peptide_metrices))
      if (n_levels == 2) {
        column_levels <- c("type", "count_level")
        leave_last_level_wide = TRUE
      } else if (!peptide_level_last) {
        column_levels <- c("type", "count_level", project_params$group_header, project_params$additional_header, replicate_info)
        leave_last_level_wide = FALSE
      } else {
        column_levels <- c("type", paste0("extra_level_", c(1:max(1, n_levels-2))), "count_level")
        leave_last_level_wide = TRUE
      }
      counts <- final_table %>%
        prepare_final_table_section(columns_to_pick = peptide_count_columns, column_levels = column_levels,
                                    leave_last_level_wide = leave_last_level_wide, drop_levels = c("type")) %>%
        dplyr::filter_at(dplyr::vars(any_of(c("value", "peptides", "msms_count", "ms_ms_count"))), ~ !is.na(.x))
      
      if (!peptide_level_last & n_levels > 2) {
        counts <- counts %>%
          tidyr::unite(group, project_params$group_header, remove = !keep_group_header) %>%
          dplyr::mutate_at(project_params$additional_header, 
                           \(x) as.numeric(gsub(project_params$additional_header_unit, "", x))) %>% 
          tidyr::pivot_wider(names_from = count_level, values_from = value) %>%
          janitor::clean_names()
      } else if (peptide_level_last & n_levels > 2) {
        counts <- counts %>%
          tidyr::unite(peptide_group, paste0("extra_level_", c(1:max(1, n_levels-2))))
      }
    counts
    },
    error = function (e) {
      if (e$message != "argument is of length zero") {
        message(paste(e, "MSMS/Peptide counts could not be processed."), sep = "\n")
      }
      NULL
    }
  )
  
  # (4) INTENSITIES
  intensity_regex <- "^(Batch-corrected|Imputed|(Unn|N)ormali[sz]ed)(.*)[Ii]ntensit(.*)"
  intensity_columns <- regmatches(final_table_columns, regexpr(intensity_regex, final_table_columns)) 
  
  intensities <- tryCatch(
    expr = {
      n_levels <- unique(purrr::map_int(intensity_columns, \(x) length(stringr::str_split_1(x, "\\|"))))
      column_levels <- c("type", project_params$group_header, project_params$additional_header, replicate_info)
      if (n_levels > length(column_levels)) {
        column_levels <- append(column_levels, extra_level_info, after = 1)
      }
      final_table %>%
        prepare_final_table_section(columns_to_pick = intensity_columns, column_levels = column_levels,
                                    leave_last_level_wide = FALSE, drop_levels = c()) %>%
        tidyr::unite(group, project_params$group_header, remove = !keep_group_header) %>%
        dplyr::mutate_at(project_params$additional_header, 
                         \(x) as.numeric(gsub(project_params$additional_header_unit, "", x))) %>%
        dplyr::mutate(type = tolower(stringr::str_replace_all(type, c("(LFQ|(TMT|tmt) reporter)" = "", 
                                                                      "\\s+" = "_", 
                                                                      "[Ii]ntensities" = "intensity"))), 
                      .after = scout_id) %>%
        tidyr::pivot_wider(names_from = type, values_from = value)
      
    },
    error = function (e) {
      if (e$message != "argument is of length zero") {
        message(paste(e, "Intensities could not be processed."), sep = "\n")
      }
      NULL
    }
  )
  
  # (5) SCALED INTENSITIES
  scaled_intensity_regex <- "^(Scale|Rel)(.*)intensit(.*)|^Intensity ratio(.*)"
  scaled_intensity_columns <- regmatches(final_table_columns, regexpr(scaled_intensity_regex, final_table_columns)) 
  
  scaled_intensities <- tryCatch(
    expr = {
      n_levels <- unique(purrr::map_int(scaled_intensity_columns, \(x) length(stringr::str_split_1(x, "\\|"))))
      column_levels <- c("type", project_params$group_header, project_params$additional_header, replicate_info)
      if (n_levels > length(column_levels)) {
        column_levels <- append(column_levels, extra_level_info, after = 1)
      }
      final_table %>%
        prepare_final_table_section(columns_to_pick = scaled_intensity_columns, column_levels = column_levels,
                                    leave_last_level_wide = FALSE, drop_levels = c("type")) %>%
        tidyr::unite(group, project_params$group_header, remove = !keep_group_header) %>%
        dplyr::mutate_at(project_params$additional_header, 
                         \(x) as.numeric(gsub(project_params$additional_header_unit, "", x))) %>%
        dplyr::rename(scaled_intensity = value)
    },
    error = function (e) {
      if (e$message != "argument is of length zero") {
        message(paste(e, "Scaled intensities could not be processed."), sep = "\n")
      }
      NULL
    }
  )
  
  # (6) CURVE PARAMETERS
  error_idx <- grep("(^Curve parameters(.*)|^Curve fitting(.*))Standard [Ee]rror", final_table_columns)
  final_table_columns[error_idx] <- paste0(final_table_columns[error_idx], "_error")
  colnames(final_table) <- final_table_columns
  
  curve_parameter_regex <- "^Curve parameters(.*)|^Curve fitting(.*)"
  curve_parameter_columns <- regmatches(final_table_columns, regexpr(curve_parameter_regex, final_table_columns))

  curve_parameters <- tryCatch(
    expr = {
      n_levels <- unique(purrr::map_int(curve_parameter_columns, \(x) length(stringr::str_split_1(x, "\\|"))))
      column_levels <- c("type", if (project_params$quant_type == "maxquant_lfq") extra_level_info else NULL, 
                         project_params$group_header, "parameter")
      if (n_levels > length(column_levels)) {
        column_levels <- append(column_levels, "parameter_group", after = length(column_levels) - 1)
      }
      final_table %>%
        prepare_final_table_section(columns_to_pick = curve_parameter_columns, column_levels = column_levels,
                                    leave_last_level_wide = TRUE, 
                                    drop_levels = intersect(column_levels, c("type", "parameter_group"))) %>%
        tidyr::unite(group, project_params$group_header, remove = !keep_group_header) %>%
        add_pec50_column(concentration_unit = project_params$additional_header_unit)
    },
    error = function (e) {
      if (e$message != "argument is of length zero") {
        message(paste(e, "Curve parameters could not be processed."), sep = "\n")
      }
      NULL
    }
  )
  
  # (7) CR PARAMETERS
  cr_parameter_regex <- "^Mean competition ratio(.*)"
  cr_parameter_columns <- regmatches(final_table_columns, regexpr(cr_parameter_regex, final_table_columns))
  
  cr_parameters <- tryCatch(
    expr = {
      column_levels <- c("type", project_params$group_header)
      final_table %>%
        prepare_final_table_section(columns_to_pick = cr_parameter_columns, column_levels = column_levels,
                                    leave_last_level_wide = TRUE, drop_levels = c("type")) %>%
        tidyr::unite(group, project_params$group_header, remove = !keep_group_header)
    },
    error = function (e) {
      if (e$message != "`cols` must select at least one column." & e$message != "") {
        message(paste(e, "CR parameters could not be processed."), sep = "\n")
      }
      NULL
    }
  )
  
  # (8) STATISTICS - T-TEST
  t_test_regex <- "^Statistical testing results(.*)Pairwise comparisons \\(t-test\\)(.*)"
  t_test_columns <- regmatches(final_table_columns, regexpr(t_test_regex, final_table_columns))

  statistics_t <- tryCatch(
    expr = {
      column_levels <- c("type", "description", "test_type", "comparison", "parameter")
      final_table %>%
        prepare_final_table_section(columns_to_pick = t_test_columns, column_levels = column_levels,
                                    leave_last_level_wide = TRUE, drop_levels = c("type", "test_type"))
    },
    error = function (e) {
      if (e$message != "`cols` must select at least one column.") {
        message(paste(e, "T-test statistics could not be processed."), sep = "\n")
      }
      NULL
    }
  )
  
  # (9) STATISTICS - F-TEST
  f_test_regex <- "^Statistical testing results(.*)Any effect \\(ANOVA\\)(.*)"
  f_test_columns <- regmatches(final_table_columns, regexpr(f_test_regex, final_table_columns))
  
  statistics_F <- tryCatch(
    expr = {
      column_levels <- c("type", "description", "test_type", "parameter")
      final_table %>%
        prepare_final_table_section(columns_to_pick = f_test_columns, column_levels = column_levels,
                                    leave_last_level_wide = TRUE) %>%
        dplyr::select(-type, -test_type)
    },
    error = function (e) {
      if (e$message != "`cols` must select at least one column.") {
        message(paste(e, "F-test statistics could not be processed."), sep = "\n")
      }
      NULL
    }
  )
  
  # (10) STATISTICS - F-TEST CURVES
  f_curve_regex <- "^F-Test results(.*)|^Statistical testing results for denaturation profiles(.*)"
  f_curve_columns <- regmatches(final_table_columns, regexpr(f_curve_regex, final_table_columns))
  
  curve_statistics_F <- tryCatch(
    expr = {
      column_levels <- c("type", "comparison", "parameter")
      final_table %>%
        prepare_final_table_section(columns_to_pick = f_curve_columns, column_levels = column_levels,
                                    leave_last_level_wide = TRUE) %>%
        dplyr::filter_at(dplyr::vars(any_of(c("f_statistic", "p_value_uncorrected"))), ~ !is.na(.x)) %>%
        dplyr::select(-type)
    },
    error = function (e) {
      if (e$message != "`cols` must select at least one column.") {
        message(paste(e, "Curve parameter F-test statistics could not be processed."), sep = "\n")
      }
      NULL
    }
  )
  
  return(list("identifier" = identifier, "annotations" = annotations, "peptide_counts" = peptide_counts,
              "intensities" = intensities, "scaled_intensities" = scaled_intensities, 
              "curve_parameters" = curve_parameters, "cr_parameters" = cr_parameters, 
              "statistics_t" = statistics_t, "statistics_F" = statistics_F, "curve_statistics_F" = curve_statistics_F))
}


#' Filter all subtables by a subset of scout IDs 
#'
#' @param data project data, as result of prepare_final_table()
#' @param scout_id_list list of scout IDs that should remain in the data set
#' @return an updated version of the input data with filtered scout_id column
#'
#' @examples
#' pois <- c("KRAS_P01116", "KRAS_P01116-2")
#' filtered_data <- filter_data_by_scout_ids(data, pois)
#' filtered_data <- data %>% filter_data_by_scout_ids(pois)
filter_data_by_scout_ids <- function(data, scout_id_list) {
  
  filtered_data <- purrr::map2(data, names(data), \(x, y) {
    tryCatch({x %>% dplyr::filter(scout_id %in% scout_id_list)},
             error = function (e) {tibble::tibble()})
  })
  
  return(filtered_data)
}


#' Add pEC50 to the curve parameters with a -log10 transform
#' 
#' If the concentration unit is uM or nM the pEC50 value is converted to M, else there is no further unit change
#'
#' @param curve_parameters curve parameters from final table parsed by prepare_final_table() function
#' @param concentration_unit concentration unit of EC50 value
#' @return updated curve parameter table with pec50 column
#'
#' @examples curve_parameters %>% add_pec50_column(concentration_unit = project_params$additional_header_unit)
add_pec50_column <- function(curve_parameters, concentration_unit) {
  
  if ("ec50" %in% names(curve_parameters)) {
    
    # depending on the unit of the EC50 the pEC50 needs to be corrected to be in M (standard unit for pEC50)
    # due to the log10 transform the correction factor is the change in the order of magnitude
    concentration_unit <- if (is.null(concentration_unit)) "-" else concentration_unit
    if (concentration_unit == "nM") {
      correction_factor <- 9
    } else if (concentration_unit == "uM") {
      correction_factor <- 6
    } else {
      correction_factor <- 0
    }
    # pEC50 [M] = -log10(EC50 [M])
    updated_curve_parameters <- curve_parameters %>%
      dplyr::mutate(pec50 = -log10(ec50) + correction_factor, .after = "ec50")
    
    message(paste("The pEC50 in", 
                  if (correction_factor == 0) paste0("[", concentration_unit, "]") else "[M]", 
                  "was added to the curve parameters.\n"))
    
  } else {
    updated_curve_parameters <- curve_parameters
  }
  
  return(updated_curve_parameters)
}
