# Degradomics reporting demo dashboard - synthetic data only.
# Files in R/ are sourced automatically by Shiny before this file runs.

library(shiny)
library(bslib)
library(magrittr)

demo_data <- make_demo_data()
colours <- condition_colours(unique(demo_data$samples$compound))
fold_changes <- sample_fold_changes(demo_data$intensities)

synthetic_banner <- div(
  class = "alert alert-warning py-2 mb-3",
  strong("Synthetic data. "),
  "All compounds, samples and measurements are simulated for demonstration; protein names other than ",
  "the public neosubstrates are made up (SYNG..., SYNZF...)."
)

ui <- page_navbar(
  title = "ProteomeScout HT - showcase",
  theme = bs_theme(version = 5, primary = "#11447b", secondary = "#76afe2", danger = "#c40a4f"),
  fillable = FALSE,
  header = synthetic_banner,
  nav_panel("Screen status", mod_screen_ui("screen")),
  nav_panel("Experimental setup", mod_setup_ui("setup")),
  nav_panel("QC", mod_qc_ui("qc")),
  nav_panel("Statistics", mod_statistics_ui("statistics")),
  nav_panel("Proteins of interest", mod_proteins_ui("proteins"))
)

server <- function(input, output, session) {
  mod_screen_server("screen", demo_data$screen)
  mod_setup_server("setup", demo_data$samples, colours)
  mod_qc_server("qc", demo_data$intensities, demo_data$samples, colours)
  regulated_stats <- mod_statistics_server(
    "statistics", demo_data$statistics, demo_data$proteins, fold_changes, demo_data$samples$sample
  )
  mod_proteins_server("proteins", demo_data$proteins, fold_changes, regulated_stats, colours)
}

shinyApp(ui, server)
