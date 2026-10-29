# Overview: entry page for the study's client - is the screen on track (progress, positive control, plates) and
# what came out of it (compound card wheel). Technical QC lives in the QC tab.

mod_overview_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h5(class = "ov-section", "Study progress", span(class = "ov-sub", "Degrader screen · 3 × 96-well plates")),
    uiOutput(ns("tiles")),
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("Runs acquired against plan", div(class = "ov-muted", "Cumulative runs, 96 per plate")),
        plotly::plotlyOutput(ns("runs"), height = "300px"),
        card_footer(class = "ov-muted", "Solid: acquired. Dashed: plan.")
      ),
      card(
        card_header(
          "Positive control recovery per plate",
          div(class = "ov-muted", "GSPT1 reduction by the positive control, both replicates: the assay works")
        ),
        plotly::plotlyOutput(ns("pos_ctrl"), height = "300px"),
        card_footer(class = "ov-muted", textOutput(ns("pos_ctrl_note"), inline = TRUE))
      )
    ),
    div(
      class = "ov-section d-flex align-items-center",
      h5(class = "m-0", "Compound highlights"),
      span(class = "ov-sub", "Compounds degrading known neosubstrates, most first"),
      div(
        class = "ms-auto d-flex align-items-center gap-2",
        tags$button(class = "btn btn-sm btn-outline-secondary", onclick = sprintf("ovWheel('%s', -1)", ns("wheel")),
                    `aria-label` = "Previous compounds", icon("chevron-left")),
        span(class = "ov-muted", id = ns("wheel_pos")),
        tags$button(class = "btn btn-sm btn-outline-secondary", onclick = sprintf("ovWheel('%s', 1)", ns("wheel")),
                    `aria-label` = "Next compounds", icon("chevron-right"))
      )
    ),
    uiOutput(ns("cards")),
    p(class = "ov-note", textOutput(ns("card_note"), inline = TRUE)),
    h5(class = "ov-section", "Plates"),
    card(tableOutput(ns("plates")))
  )
}

badge_tag <- function(badge) {
  span(class = paste("ov-badge", paste0("ov-", badge$level)), badge$label)
}

reduction_row <- function(gene, reduction, hit) {
  div(
    class = "ov-red",
    span(class = "ov-red-gene", gene),
    span(class = "ov-red-track", span(
      class = if (hit) "ov-red-bar ov-hit" else "ov-red-bar",
      style = sprintf("width: %d%%;", as.integer(reduction))
    )),
    span(class = "ov-red-value", paste0(reduction, " %"))
  )
}

compound_card_tag <- function(card, ns) {
  top_rows <- if (nrow(card$top) == 0) {
    div(class = "ov-muted", "No neosubstrate reduced")
  } else {
    purrr::pmap(card$top, \(gene_name, hit, reduction) reduction_row(gene_name, reduction, hit))
  }
  neo_label <- "known neosubstrates degraded"
  if (card$n_candidates > 0) {
    plural <- if (card$n_candidates > 1) "s" else ""
    neo_label <- sprintf("%s (+%d new candidate%s)", neo_label, card$n_candidates, plural)
  }
  div(
    class = "card ov-card",
    div(
      class = "card-body",
      div(
        class = "d-flex justify-content-between align-items-baseline",
        strong(card$compound),
        span(class = "ov-muted", sprintf("plate %s · %s", sub("P", "", card$plate_id), card$wells))
      ),
      div(lapply(card$badges, badge_tag)),
      div(
        class = "d-flex gap-4",
        div(div(class = "ov-number", card$n_neo), div(class = "ov-muted", neo_label)),
        div(div(class = "ov-number", card$n_degraded), div(class = "ov-muted", "proteins degraded in total"))
      ),
      div(class = "ov-muted", "Strongest reduction among neosubstrates"),
      top_rows,
      tags$a(
        class = "ov-link", href = "#",
        onclick = sprintf(
          "Shiny.setInputValue('%s', '%s', {priority: 'event'}); return false;", ns("open_compound"), card$compound
        ),
        "Details in Statistics", icon("arrow-right")
      )
    )
  )
}

#' @return reactive with the compound whose card link was clicked (NULL before)
mod_overview_server <- function(id, study, today = Sys.Date()) {
  moduleServer(id, function(input, output, session) {
    stats <- annotated_statistics(study)
    progress <- progress_summary(study, stats)
    plates <- plate_status(study)$plate_id

    output$tiles <- renderUI({
      tile <- function(label, value, detail) {
        card(div(class = "ov-muted", label), div(class = "ov-big", value), div(class = "ov-muted", detail))
      }
      layout_columns(
        col_widths = c(3, 3, 3, 3),
        tile(
          "Plates", sprintf("%d / %d", progress$n_plates_done, progress$n_plates),
          sprintf("%d of %d runs acquired", progress$runs_acquired, progress$runs_total)
        ),
        tile(
          "Runs passing QC", sprintf("%.1f %%", progress$qc_pass_rate),
          sprintf("%d runs flagged for review", progress$n_qc_flagged)
        ),
        tile(
          "Positive control", sprintf("%.0f %%", progress$pos_ctrl),
          sprintf("median %s reduction per plate", positive_control_readout)
        ),
        tile(
          "Neosubstrate hits", progress$n_neo_compounds,
          sprintf(
            "compounds degrading a known neosubstrate (%d of %d with results)",
            progress$n_neo_compounds, progress$n_compounds_results
          )
        )
      )
    })

    output$runs <- plotly::renderPlotly(runs_plan_plot(runs_against_plan(study$run_log), today))
    output$pos_ctrl <- plotly::renderPlotly(positive_control_plot(positive_control_replicates(study), plates))
    output$pos_ctrl_note <- renderText({
      replicates <- positive_control_replicates(study)
      below <- replicates %>% dplyr::filter(reduction < pos_ctrl_acceptance)
      paste(
        "Hollow circles: plates not yet acquired.",
        if (nrow(below) == 0) {
          "All replicates within acceptance."
        } else {
          sprintf("%d replicate(s) below acceptance, flagged for review.", nrow(below))
        }
      )
    })

    cards <- lapply(card_compounds(stats), compound_card, stats = stats, plate_assignment = study$plate_assignment)
    output$cards <- renderUI({
      tagList(
        div(id = session$ns("wheel"), class = "ov-wheel", lapply(cards, compound_card_tag, ns = session$ns)),
        tags$script(sprintf("ovWheelInit('%s', '%s');", session$ns("wheel"), session$ns("wheel_pos")))
      )
    })
    output$card_note <- renderText(paste(
      "Degraded = reduced by at least 50 % compared to the DMSO wells of the same plate (log2 fold change <= -1)",
      "with a BH-adjusted p-value <= 0.05. Bars show the reduction in protein amount; grey bars are not",
      "significant. Badges: 'broad' from 3 known neosubstrates; 'selective' when at least 25 % of the degraded",
      "proteins are neosubstrates; caution when fewer than 10 % of 10 or more degraded proteins are. Badge",
      "thresholds are demo conventions."
    ))

    output$plates <- renderTable(
      {
        plates_table(study, stats) %>%
          dplyr::transmute(
            Plate = sub("P", "", plate_id), Status = status, Week = week,
            `Pos. ctrl` = ifelse(is.na(reduction), "-", sprintf("%.0f %%", reduction)),
            `DMSO CV` = ifelse(is.na(dmso_cv), "-", sprintf("%.1f %%", dmso_cv)),
            `Runs passing QC` = ifelse(runs_checked == 0, "-", sprintf("%d / %d", runs_passed, n_runs)),
            `Neosubstrate hits` = ifelse(is.na(n_neo_compounds), "-", as.character(n_neo_compounds))
          )
      },
      width = "100%"
    )

    reactive(input$open_compound)
  })
}
