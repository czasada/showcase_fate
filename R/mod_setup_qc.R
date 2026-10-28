# Experimental setup and QC analysis of the deep-dive study (sections 1-2 of the WP1 degradomics report).

mod_setup_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = c(5, 7),
    card(
      card_header("Study design"),
      p(
        "Deep-dive study of five compounds from the screen, each compared to DMSO in a synthetic cell line ",
        "(4 replicates, 1 uM, 24 h). Proteins quantified by label-free DIA, median-normalized, missing values ",
        "imputed (left-censored), compared to DMSO with a Welch t-test and Benjamini-Hochberg correction."
      ),
      tableOutput(ns("conditions"))
    ),
    card(
      full_screen = TRUE,
      card_header("Plate layout (96 wells)"),
      plotly::plotlyOutput(ns("plate"), height = "360px")
    )
  )
}

mod_setup_server <- function(id, samples, colours) {
  moduleServer(id, function(input, output, session) {
    output$conditions <- renderTable({
      samples %>%
        dplyr::summarise(
          `Cell line` = dplyr::first(cell_line), Concentration = dplyr::first(compound_concentration),
          Time = dplyr::first(time), Replicates = dplyr::n(), Wells = paste(well, collapse = ", "),
          .by = compound
        ) %>%
        dplyr::rename(Condition = compound)
    })
    output$plate <- plotly::renderPlotly({
      wells <- expand.grid(plate_row = LETTERS[1:8], plate_col = 1:12, stringsAsFactors = FALSE) %>%
        dplyr::left_join(samples, by = c("plate_row", "plate_col")) %>%
        dplyr::mutate(
          compound = factor(compound, levels = names(colours)),
          text = ifelse(is.na(sample), paste0(plate_row, plate_col, ": empty"), paste0(well, ": ", sample))
        )
      plate_plot(wells, "compound", "Condition", 8, 12, colours = colours, marker_size = 30)
    })
  })
}

mod_qc_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    full_screen = TRUE,
    nav_panel(
      "ID counts",
      layout_columns(
        col_widths = c(6, 6),
        plotly::plotlyOutput(ns("id_plate"), height = "330px"),
        plotly::plotlyOutput(ns("id_bar"), height = "520px")
      )
    ),
    nav_panel(
      "Summed intensity",
      layout_columns(
        col_widths = c(6, 6),
        plotly::plotlyOutput(ns("intensity_plate"), height = "330px"),
        plotly::plotlyOutput(ns("intensity_bar"), height = "520px")
      )
    ),
    nav_panel(
      "Intensity distributions",
      radioButtons(
        ns("intensity_type"), NULL,
        c("Raw" = "raw_log2_intensity", "Normalized" = "normalized_log2_intensity"),
        selected = "normalized_log2_intensity", inline = TRUE
      ),
      plotly::plotlyOutput(ns("distributions"), height = "560px")
    ),
    nav_panel(
      "Reproducibility",
      p("Coefficients of variation per protein and condition, on measured normalized intensities."),
      plotly::plotlyOutput(ns("cv"), height = "420px")
    ),
    nav_panel(
      "PCA",
      uiOutput(ns("pca_note")),
      layout_columns(
        col_widths = c(6, 6),
        plotly::plotlyOutput(ns("pca_12"), height = "430px"),
        plotly::plotlyOutput(ns("pca_13"), height = "430px")
      )
    ),
    nav_panel(
      "Missing values",
      layout_columns(
        col_widths = c(6, 6),
        plotly::plotlyOutput(ns("completeness"), height = "560px"),
        plotly::plotlyOutput(ns("valid_values"), height = "400px")
      )
    )
  )
}

mod_qc_server <- function(id, intensities, samples, colours) {
  moduleServer(id, function(input, output, session) {
    sample_order <- samples$sample
    qc <- sample_qc_summary(intensities, samples) %>%
      dplyr::mutate(
        compound = factor(compound, levels = names(colours)),
        text_ids = paste0(sample, " (", well, ")<br>protein IDs: ", n_ids),
        text_intensity = paste0(sample, " (", well, ")<br>summed intensity: ", signif(summed_intensity, 3))
      )

    output$id_plate <- plotly::renderPlotly({
      plate_plot(qc %>% dplyr::mutate(text = text_ids), "n_ids", "Protein IDs", 8, 12, marker_size = 26)
    })
    output$id_bar <- plotly::renderPlotly({
      sample_bar_plot(qc %>% dplyr::mutate(text = text_ids), "n_ids", "Protein IDs", sample_order, colours)
    })
    output$intensity_plate <- plotly::renderPlotly({
      plate_plot(
        qc %>% dplyr::mutate(text = text_intensity), "summed_intensity", "Summed intensity", 8, 12,
        marker_size = 26
      )
    })
    output$intensity_bar <- plotly::renderPlotly({
      sample_bar_plot(
        qc %>% dplyr::mutate(text = text_intensity), "summed_intensity", "Summed raw intensity", sample_order, colours
      )
    })
    output$distributions <- plotly::renderPlotly({
      intensity_distribution_plot(intensities, input$intensity_type, sample_order, colours)
    })
    output$cv <- plotly::renderPlotly({
      cv_plot(cv_per_condition(intensities), names(colours), colours)
    })

    pca <- pca_scores(intensities, samples)
    output$pca_note <- renderUI({
      p(sprintf("PCA on %d proteins without missing values (normalized log2 intensities, centred).", pca$n_proteins))
    })
    output$pca_12 <- plotly::renderPlotly(pca_plot(pca, colours, "PC1", "PC2"))
    output$pca_13 <- plotly::renderPlotly(pca_plot(pca, colours, "PC1", "PC3"))
    output$completeness <- plotly::renderPlotly(completeness_plot(intensities, sample_order))
    output$valid_values <- plotly::renderPlotly(valid_values_plot(intensities, length(sample_order)))
  })
}
