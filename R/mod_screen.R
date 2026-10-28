# Screen status: where every well of the 384-well screen stands, independent of the deep-dive compounds.

screen_colour_choices <- c(
  "Stage" = "stage", "Protein IDs" = "n_ids",
  "Degraded proteins" = "n_degraded", "Neosubstrates degraded" = "n_neosubstrates_degraded"
)

mod_screen_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("value_boxes")),
    layout_columns(
      col_widths = c(8, 4),
      card(
        full_screen = TRUE,
        card_header(
          class = "d-flex justify-content-between align-items-center",
          "Plate map (384 wells)",
          selectInput(ns("colour_by"), NULL, screen_colour_choices, width = "220px")
        ),
        plotly::plotlyOutput(ns("plate"), height = "430px"),
        card_footer(
          "Column 1: DMSO. Column 24: positive control (A-H) and blanks (I-P). ",
          "Deep-dive compounds (detailed report tabs) are shown in bold in the table below."
        )
      ),
      card(
        full_screen = TRUE,
        card_header("Progress by stage"),
        plotly::plotlyOutput(ns("stages"), height = "430px")
      )
    ),
    layout_columns(
      col_widths = c(5, 7),
      card(
        full_screen = TRUE,
        card_header("Compounds with most degraded proteins (* deep dive)"),
        plotly::plotlyOutput(ns("top_compounds"), height = "420px")
      ),
      card(
        full_screen = TRUE,
        card_header("Compounds"),
        DT::DTOutput(ns("table"))
      )
    )
  )
}

mod_screen_server <- function(id, screen) {
  moduleServer(id, function(input, output, session) {
    compounds <- screen %>% dplyr::filter(well_type == "compound")

    output$value_boxes <- renderUI({
      n <- nrow(compounds)
      share <- function(stages) {
        k <- sum(compounds$stage %in% stages)
        sprintf("%d (%.0f %%)", k, k / n * 100)
      }
      analysed <- compounds %>% dplyr::filter(stage %in% c("analysed", "reported"))
      layout_columns(
        fill = FALSE,
        value_box("Compounds on plate", n, theme = "primary"),
        value_box("Acquired", share(c("acquired", "QC failed", "QC passed", "analysed", "reported"))),
        value_box("QC passed", share(c("QC passed", "analysed", "reported"))),
        value_box("Reported", share("reported")),
        value_box("QC failed", sum(compounds$stage == "QC failed"), theme = "warning"),
        value_box(
          "Compounds degrading a neosubstrate",
          sprintf("%d of %d analysed", sum(analysed$n_neosubstrates_degraded > 0), nrow(analysed))
        )
      )
    })

    output$plate <- plotly::renderPlotly({
      fill <- input$colour_by
      wells <- screen %>%
        dplyr::mutate(text = paste0(
          "<b>", well, "</b> ", compound, "<br>stage: ", stage,
          "<br>protein IDs: ", ifelse(is.na(n_ids), "-", n_ids),
          "<br>degraded: ", ifelse(is.na(n_degraded), "-", n_degraded),
          "<br>neosubstrates degraded: ", ifelse(is.na(n_neosubstrates_degraded), "-", n_neosubstrates_degraded),
          ifelse(is.na(top_neosubstrate), "", paste0("<br>top neosubstrate: ", top_neosubstrate)),
          ifelse(deep_dive, "<br><i>deep-dive compound</i>", "")
        ))
      if (fill == "stage") {
        plate_plot(wells, "stage", "Stage", 16, 24, colours = stage_colours, marker_size = 17)
      } else {
        plate_plot(wells, fill, names(screen_colour_choices)[screen_colour_choices == fill], 16, 24, marker_size = 17)
      }
    })

    output$stages <- plotly::renderPlotly({
      df <- compounds %>%
        dplyr::count(stage, .drop = FALSE) %>%
        dplyr::mutate(text = paste0(stage, ": ", n, " compounds"))
      plotly::plot_ly(
        df,
        x = ~n, y = ~stage, color = ~stage, colors = stage_colours, type = "bar", orientation = "h",
        text = ~n, textposition = "outside", hovertext = ~text, hoverinfo = "text"
      ) %>%
        finish_plot(xaxis = list(title = "Compounds"), yaxis = category_axis(rev(screen_stages)), showlegend = FALSE)
    })

    output$top_compounds <- plotly::renderPlotly({
      df <- compounds %>%
        dplyr::filter(!is.na(n_degraded)) %>%
        dplyr::slice_max(n_degraded, n = 20, with_ties = FALSE) %>%
        dplyr::mutate(
          label = ifelse(deep_dive, paste0(compound, " *"), compound),
          neosubstrates = factor(pmin(n_neosubstrates_degraded, 3), levels = 0:3, labels = c("0", "1", "2", "3+")),
          text = paste0(compound, "<br>degraded: ", n_degraded, "<br>neosubstrates: ", n_neosubstrates_degraded)
        )
      plotly::plot_ly(
        df,
        x = ~n_degraded, y = ~label, color = ~neosubstrates, type = "bar", orientation = "h",
        colors = c("0" = "#C8C8C8", "1" = os_colours[2], "2" = os_colours[1], "3+" = os_colours[9]),
        text = ~text, hoverinfo = "text", textposition = "none"
      ) %>%
        finish_plot(
          xaxis = list(title = "Degraded proteins"), yaxis = category_axis(rev(as.character(df$label))),
          legend = list(title = list(text = "Neosubstrates<br>degraded"))
        )
    })

    output$table <- DT::renderDT({
      df <- compounds %>%
        dplyr::arrange(dplyr::desc(deep_dive), dplyr::desc(n_degraded)) %>%
        dplyr::transmute(
          Compound = compound, Well = well, Stage = stage, `Protein IDs` = n_ids, `Degraded` = n_degraded,
          `Neosubstrates degraded` = n_neosubstrates_degraded, `Top neosubstrate` = top_neosubstrate,
          `Deep dive` = ifelse(deep_dive, "yes", "")
        )
      DT::datatable(
        df,
        rownames = FALSE, filter = "top", selection = "none",
        options = list(pageLength = 8, dom = "tip", scrollX = TRUE)
      ) %>%
        DT::formatStyle("Deep dive", target = "row", fontWeight = DT::styleEqual("yes", "bold"))
    })
  })
}
