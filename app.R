# Degradomics reporting demo dashboard - synthetic data only.
# Files in R/ are sourced automatically by Shiny before this file runs.

library(shiny)
library(bslib)
library(magrittr)

study <- load_demo_study()

synthetic_note <- paste(
  "Synthetic data: all compounds, samples and measurements are simulated for demonstration; protein names other",
  "than the public neosubstrates are made up (SYNG..., SYNZF...)."
)

ui <- page_navbar(
  id = "main",
  title = "Showcase FATE",
  theme = bs_theme(version = 5, primary = "#11447b", secondary = "#76afe2", danger = "#c40a4f"),
  fillable = FALSE,
  header = tags$head(
    tags$link(rel = "stylesheet", href = "overview.css"),
    tags$script(src = "overview.js")
  ),
  footer = div(class = "text-muted small text-center py-3", synthetic_note),
  nav_panel("Overview", mod_overview_ui("overview")),
  nav_panel("Experimental setup", mod_setup_ui("setup")),
  nav_panel("QC", mod_qc_ui("qc")),
  nav_panel("Statistics", mod_statistics_ui("statistics")),
  nav_panel("Proteins of interest", mod_proteins_ui("proteins")),
  nav_spacer(),
  nav_item(span(class = "badge text-bg-warning", title = synthetic_note, "Synthetic data"))
)

server <- function(input, output, session) {
  open_compound <- mod_overview_server("overview", study)
  observeEvent(open_compound(), nav_select("main", "Statistics"))
  mod_setup_server("setup", study)
  mod_qc_server("qc", study)
  mod_statistics_server("statistics", study, focus_compound = open_compound)
  mod_proteins_server("proteins", study)
}

shinyApp(ui, server)
