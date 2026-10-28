# Statistical analysis (section 3 of the WP1 degradomics report): each compound vs DMSO.
# Default cutoffs follow the report header: BH-adjusted p <= 0.05, |log2 FC| >= 1, both-sided.

mod_statistics_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      title = "Cutoffs",
      radioButtons(
        ns("p_col"), "p-value",
        c("BH-adjusted" = "p_value_corrected", "Uncorrected" = "p_value_uncorrected")
      ),
      numericInput(ns("p_cutoff"), "p-value cutoff (<=)", 0.05, min = 0.001, max = 1, step = 0.01),
      numericInput(ns("log2fc_cutoff"), "|log2 FC| cutoff (>=)", 1, min = 0, max = 5, step = 0.1),
      radioButtons(
        ns("direction"), "Direction",
        c("Both-sided" = "both-sided", "Only down (degraded)" = "only_down", "Only up" = "only_up")
      ),
      numericInput(ns("top_n"), "Labelled proteins per volcano", 5, min = 0, max = 30, step = 1),
      checkboxInput(ns("highlight_neo"), "Highlight established neosubstrates", TRUE)
    ),
    navset_card_tab(
      full_screen = TRUE,
      nav_panel(
        "Significant IDs",
        plotly::plotlyOutput(ns("sig_counts"), height = "380px")
      ),
      nav_panel(
        "Volcano plots",
        selectInput(ns("volcano_compound"), "Comparison", choices = NULL),
        plotly::plotlyOutput(ns("volcano"), height = "560px")
      ),
      nav_panel(
        "Overlap",
        layout_columns(
          col_widths = c(5, 7),
          plotly::plotlyOutput(ns("overlap"), height = "420px"),
          DT::DTOutput(ns("intersections"))
        )
      ),
      nav_panel(
        "Clustermap",
        p("Per-sample log2 fold changes to mean DMSO of the significant proteins with the lowest p-values (max. 100), ",
          "rows clustered (complete linkage, Euclidean distance)."),
        plotly::plotlyOutput(ns("clustermap"), height = "900px")
      ),
      nav_panel(
        "Results table",
        DT::DTOutput(ns("results"))
      )
    )
  )
}

#' @return reactive with the statistics including the regulation call under the current cutoffs
mod_statistics_server <- function(id, statistics, proteins, fold_changes, sample_order) {
  moduleServer(id, function(input, output, session) {
    stats_annotated <- statistics %>%
      dplyr::left_join(proteins %>% dplyr::select(protein_id, gene_name, neosubstrate_category), by = "protein_id")
    compounds <- unique(statistics$compound)
    updateSelectInput(session, "volcano_compound", choices = compounds)

    regulated_stats <- reactive({
      req(input$p_cutoff, !is.na(input$log2fc_cutoff))
      add_regulation(stats_annotated, input$p_col, input$p_cutoff, input$log2fc_cutoff, input$direction)
    })
    hits <- reactive(regulated_stats() %>% dplyr::filter(regulation != "none"))

    output$sig_counts <- plotly::renderPlotly({
      counts <- hits() %>%
        dplyr::count(compound, regulation) %>%
        tidyr::complete(compound = compounds, regulation = c("down", "up"), fill = list(n = 0))
      significant_count_plot(counts, compounds)
    })

    output$volcano <- plotly::renderPlotly({
      req(input$volcano_compound)
      all <- regulated_stats()
      # common axes across comparisons for comparability, as in the report
      x_max <- max(abs(all$log2_fold_change)) * 1.05
      y_max <- max(-log10(all[[input$p_col]])) * 1.05
      neosubstrates <- proteins$gene_name[proteins$neosubstrate_category == "established"]
      highlight <- if (isTRUE(input$highlight_neo)) neosubstrates else character()
      volcano_plot(
        all %>% dplyr::filter(compound == input$volcano_compound),
        input$p_col, input$p_cutoff, input$log2fc_cutoff, input$direction, input$top_n, highlight,
        c(-x_max, x_max), y_max
      )
    })

    output$overlap <- plotly::renderPlotly({
      validate(need(nrow(hits()) > 0, "No significant proteins under the current cutoffs."))
      overlap_heatmap(pairwise_overlap(hits(), compounds), compounds)
    })
    output$intersections <- DT::renderDT({
      df <- hit_intersections(hits(), proteins) %>%
        dplyr::arrange(dplyr::desc(n_compounds), compounds) %>%
        dplyr::transmute(
          Protein = gene_name, Regulation = regulation, Compounds = compounds, `n compounds` = n_compounds,
          Neosubstrate = neosubstrate_category
        )
      DT::datatable(
        df,
        rownames = FALSE, filter = "top", selection = "none", options = list(pageLength = 10, dom = "tip")
      )
    })

    output$clustermap <- plotly::renderPlotly({
      ids <- hits() %>%
        dplyr::slice_min(.data[[input$p_col]], n = 1, with_ties = FALSE, by = protein_id) %>%
        dplyr::slice_min(.data[[input$p_col]], n = 100, with_ties = FALSE) %>%
        dplyr::pull(protein_id)
      validate(need(length(ids) >= 2, "Fewer than two significant proteins under the current cutoffs."))
      clustermap_plot(fold_changes, ids, proteins, sample_order)
    })

    output$results <- DT::renderDT({
      df <- regulated_stats() %>%
        dplyr::filter(regulation != "none" | neosubstrate_category != "none") %>%
        dplyr::arrange(compound, .data[[input$p_col]]) %>%
        dplyr::transmute(
          Protein = gene_name, ID = protein_id, Comparison = comparison, `log2 FC` = round(log2_fold_change, 2),
          p = signif(p_value_uncorrected, 3), `BH-adj. p` = signif(p_value_corrected, 3),
          Regulation = regulation, `Imputed values` = n_imputed, Neosubstrate = neosubstrate_category
        )
      DT::datatable(
        df,
        rownames = FALSE, filter = "top", selection = "none",
        caption = "Significant proteins under the current cutoffs plus all annotated neosubstrates.",
        options = list(pageLength = 15, dom = "tip", scrollX = TRUE)
      )
    })

    regulated_stats
  })
}
