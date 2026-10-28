# Team colours, ported from the ProteomeScout report templates (qmd-reports/R/os_theme.R).

os_colours <- c(
  "#11447b", "#76afe2", "#c40a4f", "#a6a6a6", "#6C7A05",
  "#B0C712", "#7A5411", "#C77E00", "#7A1D63", "#0864C7",
  "#C76730", "#7A3B17", "#C7A708", "#7A680B"
)
os_colours_selected <- c(
  "#11447b", "#76afe2", "#c40a4f", "#6C7A05", "#7A1D63",
  "#B0C712", "#7A5411", "#C77E00", "#0864C7",
  "#C76730", "#7A3B17", "#C7A708", "#7A680B"
)
vehicle_colour <- "#C8C8C8"
volcano_colours <- c(down = "#c40a4f", none = "#C8C8C8", up = "#11447b")

stage_colours <- c(
  planned = "#E6E6E6", acquired = "#76afe2", `QC failed` = "#C77E00",
  `QC passed` = "#0864C7", analysed = "#11447b", reported = "#6C7A05"
)

condition_colours <- function(conditions) {
  others <- setdiff(conditions, vehicle_name)
  c(stats::setNames(vehicle_colour, vehicle_name), stats::setNames(os_colours_selected[seq_along(others)], others))
}
