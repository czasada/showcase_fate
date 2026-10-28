# Proteins of interest (barplots section of the WP1 degradomics report): fold changes and
# degradation efficiencies per sample, imputed values marked.

mod_proteins_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      title = "Proteins",
      radioButtons(
        ns("group"), "Preset",
        c("Established neosubstrates" = "neo", "Regulated in any comparison" = "regulated", "Custom" = "custom")
      ),
      selectizeInput(
        ns("proteins"), "Proteins (max. 12)",
        choices = NULL, multiple = TRUE, options = list(maxItems = 12)
      ),
      radioButtons(
        ns("value"), "Show",
        c("Log2 fold change" = "log2_fold_change", "Degradation efficiency" = "degradation_efficiency_percent")
      ),
      helpText("Regulated = significant under the cutoffs set on the Statistics tab. Open points are imputed values.")
    ),
    card(
      full_screen = TRUE,
      card_header(textOutput(ns("title"))),
      uiOutput(ns("bars_ui"))
    )
  )
}

mod_proteins_server <- function(id, proteins, fold_changes, regulated_stats, colours) {
  moduleServer(id, function(input, output, session) {
    gene_choices <- proteins %>%
      dplyr::filter(protein_id %in% fold_changes$protein_id) %>%
      dplyr::arrange(neosubstrate_category != "established", gene_name) %>%
      dplyr::pull(gene_name)
    updateSelectizeInput(session, "proteins", choices = gene_choices, server = TRUE)

    preset_genes <- reactive({
      switch(input$group,
        neo = proteins$gene_name[proteins$neosubstrate_category == "established"],
        regulated = regulated_stats() %>%
          dplyr::filter(regulation != "none") %>%
          dplyr::slice_min(p_value_uncorrected, n = 1, with_ties = FALSE, by = protein_id) %>%
          dplyr::arrange(p_value_uncorrected) %>%
          dplyr::pull(gene_name),
        custom = NULL
      )
    })
    observeEvent(preset_genes(), {
      if (input$group != "custom") {
        updateSelectizeInput(
          session, "proteins",
          choices = gene_choices, selected = utils::head(preset_genes(), 12), server = TRUE
        )
      }
    })

    output$title <- renderText({
      if (input$value == "log2_fold_change") "Log2 fold change to mean DMSO" else "Degradation efficiency [%]"
    })

    output$bars_ui <- renderUI({
      n_rows <- max(1, ceiling(length(input$proteins) / 3))
      plotly::plotlyOutput(session$ns("bars"), height = paste0(260 * n_rows + 80, "px"))
    })

    output$bars <- plotly::renderPlotly({
      validate(need(length(input$proteins) > 0, "Select at least one protein."))
      df <- fold_changes %>%
        dplyr::left_join(dplyr::select(proteins, protein_id, gene_name), by = "protein_id") %>%
        dplyr::filter(gene_name %in% input$proteins) %>%
        dplyr::mutate(gene_name = factor(gene_name, levels = input$proteins))
      if (input$value == "log2_fold_change") {
        poi_bar_plot(df, "log2_fold_change", "Log2 FC", colours, names(colours))
      } else {
        poi_bar_plot(
          df, "degradation_efficiency_percent", "Degradation efficiency [%]", colours, names(colours), c(0, 100)
        )
      }
    })
  })
}
