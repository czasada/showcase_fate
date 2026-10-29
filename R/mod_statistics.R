# Statistical analysis (section 3 of the WP1 degradomics report): each compound vs the DMSO wells of its plate.
# Default cutoffs match the Overview: BH-adjusted p <= 0.05, log2 FC <= -1 (only down).

mod_statistics_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      title = "Cutoffs",
      radioButtons(
        ns("p_col"), "p-value",
        c("BH-adjusted" = "p_value_corrected", "Uncorrected" = "p_value_uncorrected")
      ),
      numericInput(ns("p_cutoff"), "p-value cutoff (<=)", hit_p_cutoff, min = 0.001, max = 1, step = 0.01),
      numericInput(ns("log2fc_cutoff"), "|log2 FC| cutoff (>=)", abs(hit_log2fc_cutoff), min = 0, max = 5, step = 0.1),
      radioButtons(
        ns("direction"), "Direction",
        c("Only down (degraded)" = "only_down", "Both-sided" = "both-sided", "Only up" = "only_up")
      ),
      numericInput(ns("top_n"), "Labelled proteins per volcano", 5, min = 0, max = 30, step = 1),
      checkboxInput(ns("highlight_neo"), "Highlight known neosubstrates", TRUE)
    ),
    navset_card_tab(
      id = ns("tabs"),
      full_screen = TRUE,
      nav_panel(
        "Summary",
        p(
          class = "ov-muted",
          "Hits per compound under the current cutoffs, most first. Click a bar or a compound name to open its volcano."
        ),
        plotly::plotlyOutput(ns("hit_counts"), height = "380px"),
        DT::DTOutput(ns("hit_table"))
      ),
      nav_panel(
        "Volcano",
        div(
          class = "d-flex flex-wrap gap-3 align-items-end",
          selectizeInput(ns("volcano_compound"), "Compound", choices = NULL, width = "260px"),
          radioButtons(
            ns("volcano_y"), "y axis", c("p-value" = "p", "Standard error" = "se"),
            inline = TRUE
          )
        ),
        layout_columns(
          col_widths = c(7, 5),
          div(
            plotly::plotlyOutput(ns("volcano"), height = "520px"),
            p(class = "ov-muted", textOutput(ns("volcano_note"), inline = TRUE))
          ),
          div(
            div(
              class = "d-flex flex-wrap justify-content-between align-items-center gap-2",
              strong("Selected protein"),
              input_switch(ns("magic"), "Add the top 3 compounds degrading it")
            ),
            plotly::plotlyOutput(ns("protein"), height = "420px"),
            p(class = "ov-muted", textOutput(ns("protein_note"), inline = TRUE))
          )
        )
      ),
      nav_panel(
        "Results table",
        DT::DTOutput(ns("results"))
      )
    )
  )
}

#' @param focus_compound reactive with a compound to open in the Volcano sub-tab (e.g. from an Overview card)
mod_statistics_server <- function(id, study, focus_compound = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    stats <- annotated_statistics(study)
    neosubstrates <- unique(stats$gene_name[stats$neosubstrate_category == "established"])

    reg <- reactive({
      req(input$p_col, input$p_cutoff, !is.na(input$log2fc_cutoff), input$direction)
      regulated_statistics(stats, input$p_col, input$p_cutoff, input$log2fc_cutoff, input$direction)
    })
    counts <- reactive(compound_hit_counts(reg()))

    # compounds ordered by hits; keep the selection when the cutoffs change
    observe({
      choices <- counts()$compound
      selected <- isolate(input$volcano_compound)
      if (is.null(selected) || !selected %in% choices) selected <- choices[1]
      updateSelectizeInput(session, "volcano_compound", choices = choices, selected = selected)
    })
    open_volcano <- function(compound) {
      updateSelectizeInput(session, "volcano_compound", selected = compound)
      nav_select("tabs", "Volcano")
    }
    observeEvent(focus_compound(), open_volcano(focus_compound()))
    # clicks are read only once their plot has been drawn (plotly warns about unregistered sources otherwise, e.g.
    # while the Volcano sub-tab has not been opened yet)
    drawn <- reactiveValues(hit_counts = FALSE, volcano = FALSE)
    clicked <- function(name) {
      reactive({
        req(drawn[[name]])
        plotly::event_data("plotly_click", source = session$ns(name))
      })
    }
    bar_click <- clicked("hit_counts")
    observeEvent(bar_click(), open_volcano(bar_click()$customdata))

    output$hit_counts <- plotly::renderPlotly({
      drawn$hit_counts <- TRUE
      hit_count_plot(counts(), session$ns("hit_counts"))
    })
    observeEvent(input$table_compound, open_volcano(input$table_compound))
    output$hit_table <- DT::renderDT({
      df <- counts() %>%
        dplyr::transmute(
          Compound = sprintf(
            "<a href='#' onclick=\"Shiny.setInputValue('%s', '%s', {priority: 'event'}); return false;\">%s</a>",
            session$ns("table_compound"), compound, compound
          ),
          Plate = sub("P", "", plate_id), `Tested proteins` = tested, Down = down, Up = up,
          `Known neosubstrates down` = neo_down,
          `Hits with imputed values` = ifelse(is.na(imputed_share), "-", sprintf("%.0f %%", imputed_share))
        )
      DT::datatable(
        df,
        rownames = FALSE, filter = "top", selection = "none", escape = -1,
        options = list(pageLength = 10, dom = "tip")
      )
    })

    reg_compound <- reactive({
      req(input$volcano_compound)
      reg() %>% dplyr::filter(compound == input$volcano_compound)
    })
    # protein shown next to the volcano: clicked point, else the compound's strongest hit
    selected_protein <- reactiveVal(NULL)
    observeEvent(input$volcano_compound, selected_protein(NULL))
    volcano_click <- clicked("volcano")
    observeEvent(volcano_click(), selected_protein(volcano_click()$customdata))
    protein_id <- reactive({
      id <- selected_protein()
      if (is.null(id) || !id %in% reg_compound()$scout_id) default_protein(reg_compound(), input$p_col) else id
    })

    output$volcano <- plotly::renderPlotly({
      drawn$volcano <- TRUE
      highlight <- if (isTRUE(input$highlight_neo)) neosubstrates else character()
      volcano_plot(
        reg_compound(), input$volcano_y, input$p_col, input$p_cutoff, input$log2fc_cutoff, input$direction,
        input$top_n, highlight,
        source = session$ns("volcano"), selected_id = protein_id()
      )
    })
    output$volcano_note <- renderText({
      if (identical(input$volcano_y, "se")) {
        paste(
          "y axis: -log10 standard error of the log2 fold change (|log2 FC / t|). High points are measured",
          "precisely; low points rest on noisy replicates. Colours follow the p-value cutoffs. Click a point."
        )
      } else {
        "Click a point to show its replicate values."
      }
    })

    shown_compounds <- reactive({
      extra <- if (isTRUE(input$magic)) top_degraders(reg(), protein_id(), input$volcano_compound) else character()
      c(input$volcano_compound, extra)
    })
    output$protein <- plotly::renderPlotly({
      req(protein_id())
      gene <- reg_compound()$gene_name[reg_compound()$scout_id == protein_id()]
      protein_replicate_plot(protein_replicates(study, protein_id(), shown_compounds()), gene)
    })
    output$protein_note <- renderText({
      row <- reg_compound() %>% dplyr::filter(scout_id == protein_id())
      req(nrow(row) == 1)
      base <- sprintf(
        "%s: log2 FC %.2f (%.0f %% reduction), BH-adj. p %s. Hollow: imputed. Bars: means.",
        input$volcano_compound, row$log2_fold_change, fold_change_to_degradation_efficiency(2^row$log2_fold_change),
        signif(row$p_value_corrected, 2)
      )
      if (isTRUE(input$magic) && length(shown_compounds()) == 1) {
        paste(base, "No other compound degrades", row$gene_name, "under the current cutoffs.")
      } else {
        base
      }
    })

    output$results <- DT::renderDT({
      df <- reg() %>%
        dplyr::filter(regulation != "none" | neosubstrate_category != "none") %>%
        dplyr::arrange(compound, .data[[input$p_col]]) %>%
        dplyr::transmute(
          Protein = gene_name, ID = scout_id, Compound = compound, Plate = sub("P", "", plate_id),
          `log2 FC` = round(log2_fold_change, 2), `Reduction` = sprintf(
            "%.0f %%", fold_change_to_degradation_efficiency(2^log2_fold_change)
          ),
          p = signif(p_value_uncorrected, 3), `BH-adj. p` = signif(p_value_corrected, 3),
          SE = signif(standard_error, 3), Regulation = regulation, `Valid values` = valid_value_string,
          Neosubstrate = neosubstrate_category
        )
      DT::datatable(
        df,
        rownames = FALSE, filter = "top", selection = "none",
        caption = "Hits under the current cutoffs plus all annotated neosubstrates and candidates.",
        options = list(pageLength = 15, dom = "tip", scrollX = TRUE)
      )
    })
  })
}
