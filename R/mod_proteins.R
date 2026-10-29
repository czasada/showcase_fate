# Proteins of interest: screen hits rolled up per protein (left) and the selected protein across all compounds
# (right) as bubble, volcano or ranked view. Replaces the report's bar plots of selected proteins.

mod_proteins_ui <- function(id) {
  ns <- NS(id)
  layout_columns(
    col_widths = c(5, 7),
    card(
      card_header("Proteins degraded by the screen", div(class = "ov-muted", textOutput(ns("summary"), inline = TRUE))),
      checkboxGroupInput(
        ns("categories"), NULL, unname(protein_categories),
        selected = unname(protein_categories), inline = TRUE
      ),
      DT::DTOutput(ns("targets")),
      card_footer(
        class = "ov-muted",
        "Compounds = number of compounds degrading the protein (BH-adjusted p <= 0.05, at least 50 % reduction).",
        "Best = strongest reduction. Category from the final table annotation."
      )
    ),
    card(
      full_screen = TRUE,
      uiOutput(ns("header")),
      div(
        class = "d-flex flex-wrap gap-3",
        radioButtons(ns("view"), NULL, c("Bubble" = "bubble", "Volcano" = "volcano", "Ranked" = "ranked"),
                     inline = TRUE),
        radioButtons(ns("metric"), NULL, c("log2 FC" = "log2fc", "Reduction %" = "reduction",
                                           "Remaining %" = "remaining"), inline = TRUE)
      ),
      p(class = "ov-muted", textOutput(ns("view_note"), inline = TRUE)),
      plotly::plotlyOutput(ns("profile"), height = "480px")
    )
  )
}

mod_proteins_server <- function(id, study) {
  moduleServer(id, function(input, output, session) {
    stats <- annotated_statistics(study)
    targets <- protein_targets(stats, study$identifier)
    abundance <- abundance_rank(study)
    annotation <- study$identifier %>%
      dplyr::select(scout_id, gene_name, uniprot_id) %>%
      dplyr::left_join(study$annotations, by = "scout_id")

    output$summary <- renderText({
      n_cat <- table(targets$category)
      paste(
        sprintf(
          "%d proteins degraded by at least one compound (%d known neosubstrates, %d candidates);",
          nrow(targets), n_cat[["known neosubstrate"]], n_cat[["candidate"]]
        ),
        "known first, then by number of compounds."
      )
    })

    shown <- reactive(targets %>% dplyr::filter(category %in% input$categories))
    output$targets <- DT::renderDT({
      df <- shown() %>%
        dplyr::transmute(
          Protein = gene_name, UniProt = uniprot_id, Category = category,
          Best = sprintf("%.0f %%", best_reduction), Compounds = n_compounds, `Best compound` = best_compound
        )
      DT::datatable(
        df,
        rownames = FALSE, filter = "none", selection = list(mode = "single", selected = 1, target = "row"),
        options = list(pageLength = 12, dom = "ftip", scrollX = TRUE, language = list(search = "Search gene or ID:"))
      )
    })

    protein_id <- reactive({
      rows <- input$targets_rows_selected
      table <- shown()
      req(nrow(table) > 0)
      if (length(rows) == 0) table$scout_id[1] else table$scout_id[rows[1]]
    })
    profile <- reactive(protein_profile(stats, protein_id()))
    mdc <- reactive(min_detectable_change(study, protein_id()))

    output$header <- renderUI({
      info <- annotation %>% dplyr::filter(scout_id == protein_id())
      rank <- abundance$top_percent[abundance$scout_id == protein_id()]
      div(
        h5(class = "mb-0", info$gene_name, span(class = "ov-muted", info$uniprot_id)),
        div(
          class = "ov-muted",
          sprintf(
            "%s · abundance rank: top %s %% · degraded by %d of %d compounds",
            protein_categories[[info$neosubstrate_category]], if (length(rank) == 1) rank else "-",
            sum(profile()$hit), nrow(profile())
          )
        )
      )
    })
    output$view_note <- renderText({
      switch(input$view,
        bubble = paste(
          "One bubble per compound in ID order (plates shaded); red = hit; bubble size = how many proteins the",
          "compound degrades in total, so small red bubbles are selective degraders. Min. detectable change =",
          "2 x SD of this protein over the DMSO wells."
        ),
        volcano = "This protein in every compound: one point per compound, the strongest hits labelled.",
        ranked = "The 25 compounds that reduce this protein most; red = hit."
      )
    })
    output$profile <- plotly::renderPlotly({
      req(input$view, input$metric)
      switch(input$view,
        bubble = protein_bubble_plot(profile(), input$metric, mdc()),
        volcano = protein_volcano_plot(profile()),
        ranked = protein_ranked_plot(profile(), input$metric, mdc())
      )
    })
  })
}
