# Experimental setup (section 1 of the WP1 degradomics report): workflow, study design, plate layout per plate and
# the sample table of all plates.

workflow_steps <- list(
  list(icon = "flask", title = "Treat cells", detail = "compound or DMSO per well"),
  list(icon = "vial", title = "Lyse and digest", detail = "96-well plate"),
  list(icon = "wave-square", title = "DIA mass spec", detail = "one run per well"),
  list(icon = "table", title = "Pipeline", detail = "final table per plate"),
  list(icon = "chart-column", title = "Dashboard", detail = "this page")
)

setup_colour_choices <- c(
  "Well type" = "well_type", "Run status" = "run_status", "Measurement order" = "measurement_order"
)

mod_setup_ui <- function(id) {
  ns <- NS(id)
  steps <- lapply(workflow_steps, \(s) {
    div(class = "su-step", icon(s$icon), div(s$title), div(class = "ov-muted", s$detail))
  })
  tagList(
    card(
      card_header("Workflow"),
      div(class = "su-steps", Reduce(\(a, b) tagList(a, icon("arrow-right", class = "su-arrow"), b), steps))
    ),
    layout_columns(
      col_widths = c(4, 8),
      card(card_header("Study design"), tableOutput(ns("design"))),
      card(
        full_screen = TRUE,
        card_header(
          class = "d-flex flex-wrap gap-2 justify-content-between align-items-center",
          "Plate layout",
          div(
            class = "d-flex flex-wrap gap-2",
            radioButtons(ns("plate"), NULL, choices = "P1", inline = TRUE),
            selectInput(ns("colour_by"), NULL, setup_colour_choices, width = "190px")
          )
        ),
        plotly::plotlyOutput(ns("plate_map"), height = "380px"),
        card_footer(class = "ov-muted", "Compounds sit in duplicate in row pairs (framed). Hover a well for details.")
      )
    ),
    card(card_header("Samples"), DT::DTOutput(ns("samples")))
  )
}

mod_setup_server <- function(id, study) {
  moduleServer(id, function(input, output, session) {
    wells <- setup_wells(study)
    plates <- unique(wells$plate_id)
    updateRadioButtons(
      session, "plate",
      choiceNames = paste("Plate", sub("P", "", plates)), choiceValues = plates, selected = plates[1], inline = TRUE
    )

    output$design <- renderTable(design_summary(study), colnames = FALSE, width = "100%")

    output$plate_map <- plotly::renderPlotly({
      req(input$plate, input$colour_by)
      plate_wells <- wells %>%
        dplyr::filter(plate_id == input$plate) %>%
        dplyr::mutate(
          text = paste0(
            "<b>", position, "</b> ", compound, " ", replicate,
            "<br>", run_status,
            ifelse(is.na(acquired_date), "", paste0("<br>acquired ", format(acquired_date, "%d %b"))),
            ifelse(is.na(measurement_order), "", paste0("<br>measurement order ", measurement_order))
          ),
          well_type = factor(well_type, levels = names(well_type_colours)),
          run_status = factor(run_status, levels = names(run_status_colours))
        )
      colours <- switch(input$colour_by,
        well_type = well_type_colours,
        run_status = run_status_colours
      )
      label <- names(setup_colour_choices)[setup_colour_choices == input$colour_by]
      p <- plate_plot(plate_wells, input$colour_by, label, 8, 12, colours = colours, marker_size = 28)
      p %>% plotly::layout(shapes = replicate_pair_shapes(plate_wells, rev(LETTERS[1:8])))
    })

    output$samples <- DT::renderDT({
      df <- wells %>%
        dplyr::arrange(plate_id, plate_col, plate_row) %>%
        dplyr::transmute(
          Plate = sub("P", "", plate_id), Sample = sample, Compound = compound, Replicate = replicate, Well = position,
          Order = measurement_order, Planned = format(planned_date, "%d %b"),
          Acquired = ifelse(is.na(acquired_date), "", format(acquired_date, "%d %b")), QC = run_status
        )
      DT::datatable(
        df,
        rownames = FALSE, filter = "top", selection = "none",
        options = list(pageLength = 12, dom = "tip", scrollX = TRUE)
      )
    })
  })
}
