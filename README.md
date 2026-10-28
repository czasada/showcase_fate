# Degradomics screen demo app

Shiny + bslib dashboard on synthetic data (generated at start-up, seed 640, ~0.5 s).
Everything this app needs is in this folder, so it can be copied to its own GitHub repository.

## Run locally

```r
shiny::runApp("app")
```

Packages: shiny, bslib, plotly, DT, dplyr, tidyr, tibble, purrr, magrittr.
Tests (from the project folder): `testthat::test_dir("tests/testthat")`.

## Deploy to Posit Connect Cloud

Connect Cloud publishes from a GitHub repository and needs a `manifest.json` for R content
(see https://docs.posit.co/connect-cloud/ - check the current steps there).

1. In an R session that has the packages above installed: `rsconnect::writeManifest(appDir = "app")`.
2. Copy the contents of `app/` (including `manifest.json`) to the GitHub repository.
3. In Connect Cloud: Publish > Shiny > pick the repository and branch, primary file `app.R`.

## Layout

| File | Content |
|---|---|
| `app.R` | UI (page_navbar) and server wiring |
| `R/fcn_synthetic_data.R` | proteins, samples, intensities, statistics, 384-well screen |
| `R/fcn_analysis.R` | regulation calls, fold changes, QC summaries, PCA, overlaps |
| `R/fcn_plots.R` | plotly plot functions |
| `R/mod_*.R` | one module per tab |
| `R/os_colours.R` | team colours |
