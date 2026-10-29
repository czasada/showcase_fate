# QC analysis (section 2 of the WP1 degradomics report) for all plates: run quality, intensities, reproducibility
# and PCA. Technical page; the Overview carries the client-facing summary.

qc_metric_choices <- c("Protein IDs" = "n_ids", "Summed raw intensity" = "summed_intensity")

plate_choices <- function(plates) stats::setNames(plates, paste("Plate", sub("P", "", plates)))

mod_qc_ui <- function(id) {
  ns <- NS(id)
  navset_card_tab(
    full_screen = TRUE,
    nav_panel(
      "Run quality",
      div(
        class = "d-flex flex-wrap gap-3 align-items-center",
        radioButtons(ns("metric"), NULL, qc_metric_choices, inline = TRUE),
        span(class = "ov-muted", textOutput(ns("rule"), inline = TRUE))
      ),
      plotly::plotlyOutput(ns("run_strip"), height = "320px"),
      layout_columns(
        col_widths = c(7, 5),
        div(
          selectInput(ns("map_plate"), NULL, choices = NULL, width = "160px"),
          plotly::plotlyOutput(ns("run_map"), height = "330px")
        ),
        div(strong("Flagged runs"), tableOutput(ns("flagged")), tableOutput(ns("plate_summary")))
      )
    ),
    nav_panel(
      "Intensities",
      div(
        class = "d-flex flex-wrap gap-3",
        selectInput(ns("box_plate"), NULL, choices = NULL, width = "160px"),
        radioButtons(
          ns("intensity_type"), NULL, c("Raw" = "unnormalized_intensity", "Normalized" = "normalized_intensity"),
          selected = "normalized_intensity", inline = TRUE
        )
      ),
      plotly::plotlyOutput(ns("boxes"), height = "380px"),
      p(class = "ov-muted mt-2", "Normalization factor per run (column-wise median normalization), all plates:"),
      plotly::plotlyOutput(ns("norm_factors"), height = "260px")
    ),
    nav_panel(
      "Reproducibility",
      p(class = "ov-muted", "CV per protein of the six DMSO wells of each plate (proteins with at least 3 values)."),
      plotly::plotlyOutput(ns("cv"), height = "400px"),
      tableOutput(ns("cv_table"))
    ),
    nav_panel(
      "PCA",
      div(
        class = "d-flex flex-wrap gap-3 align-items-center",
        radioButtons(ns("pca_colour"), NULL, c("Plate" = "plate_id", "Well type" = "well_type"), inline = TRUE),
        span(class = "ov-muted", textOutput(ns("pca_note"), inline = TRUE))
      ),
      layout_columns(
        col_widths = c(5, 5, 2),
        plotly::plotlyOutput(ns("pca_12"), height = "380px"),
        plotly::plotlyOutput(ns("pca_13"), height = "380px"),
        plotly::plotlyOutput(ns("scree"), height = "380px")
      )
    )
    # "Missing values" sub-tab hidden until it is reworked (completeness() and valid_values_plot() are kept)
  )
}

mod_qc_server <- function(id, study) {
  moduleServer(id, function(input, output, session) {
    runs <- run_qc(study)
    plates <- plate_choices(sort(unique(runs$plate_id)))
    updateSelectInput(session, "map_plate", choices = plates)
    updateSelectInput(session, "box_plate", choices = plates)

    output$rule <- renderText(paste(
      sprintf("QC rule: a run is flagged when its protein IDs are below %.0f %%", qc_min_share * 100),
      "of the plate median (dashed). Flagged runs stay in the final table."
    ))
    output$run_strip <- plotly::renderPlotly({
      req(input$metric)
      label <- names(qc_metric_choices)[qc_metric_choices == input$metric]
      threshold <- if (input$metric == "n_ids") qc_min_share else NULL
      run_metric_plot(runs, input$metric, label, threshold)
    })
    output$run_map <- plotly::renderPlotly({
      req(input$map_plate, input$metric)
      label <- names(qc_metric_choices)[qc_metric_choices == input$metric]
      wells <- runs %>%
        dplyr::filter(plate_id == input$map_plate) %>%
        dplyr::mutate(
          plate_row = substr(position, 1, 1), plate_col = as.integer(substr(position, 2, 3)),
          text = sprintf(
            "<b>%s</b> %s<br>%s: %s<br>QC: %s", position, sample, label,
            format(signif(.data[[input$metric]], 4), big.mark = ","), qc_status
          )
        )
      plate_plot(wells, input$metric, label, 8, 12, marker_size = 26)
    })
    output$flagged <- renderTable({
      runs %>%
        dplyr::filter(qc_status == "flagged") %>%
        dplyr::transmute(
          Plate = sub("P", "", plate_id), Well = position, Sample = sample, IDs = as.integer(n_ids),
          `% of median` = sprintf("%.0f %%", share_of_median * 100)
        )
    })
    output$plate_summary <- renderTable({
      runs %>%
        dplyr::summarise(
          `Median IDs` = as.integer(stats::median(n_ids)),
          `Runs passing` = sprintf("%d / %d", sum(qc_status == "passed"), dplyr::n()),
          .by = plate_id
        ) %>%
        dplyr::mutate(Plate = sub("P", "", plate_id), .before = 1) %>%
        dplyr::select(-plate_id)
    })

    output$boxes <- plotly::renderPlotly({
      req(input$box_plate, input$intensity_type)
      label <- if (input$intensity_type == "normalized_intensity") "Normalized log2 intensity" else "Raw log2 intensity"
      intensity_box_plot(intensity_box_stats(study, input$box_plate, input$intensity_type), label)
    })
    output$norm_factors <- plotly::renderPlotly(normalization_factor_plot(normalization_factors(study)))

    cvs <- dmso_cv(study)
    output$cv <- plotly::renderPlotly(dmso_cv_plot(cvs))
    output$cv_table <- renderTable({
      cvs %>%
        dplyr::summarise(median = stats::median(cv), n = dplyr::n(), .by = c(plate_id, type)) %>%
        tidyr::pivot_wider(names_from = type, values_from = c(median, n)) %>%
        dplyr::transmute(
          Plate = sub("P", "", plate_id), `Median CV raw` = sprintf("%.1f %%", median_raw),
          `Median CV normalized` = sprintf("%.1f %%", median_normalized), Proteins = n_normalized
        )
    })

    pca <- qc_pca(study)
    output$pca_note <- renderText(sprintf(
      "PCA on %d proteins measured in every run (normalized log2 intensities, centred).", pca$n_proteins
    ))
    output$pca_12 <- plotly::renderPlotly(qc_pca_plot(pca, input$pca_colour, "PC1", "PC2"))
    output$pca_13 <- plotly::renderPlotly(qc_pca_plot(pca, input$pca_colour, "PC1", "PC3"))
    output$scree <- plotly::renderPlotly(scree_plot(pca$variance))
  })
}
