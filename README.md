# Degradomics screen demo app

Shiny + bslib dashboard on a synthetic 3 × 96-well degrader screen (seed 640). The app runs no analysis: at
start-up it writes synthetic pipeline output to a temporary folder (one flat-header final table and sample annotation
per measured plate, plus a run log) and reads it back with the reports' `prepare_final_table()` and
`infer_plate_layout()`, as it would read real pipeline output. Start-up takes about 10-15 s.
Everything this app needs is in this folder, so it can be copied to its own GitHub repository.

## Run locally

```r
shiny::runApp("app")
```

Packages: shiny, bslib, plotly, DT, dplyr, tidyr, tibble, purrr, magrittr, stringr, data.table, janitor.
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
| `R/fcn_synthetic_data.R` | synthetic study: per-plate final tables, sample annotations, run log |
| `R/fcn_final_table.R` | final table reader, copied unchanged from the qmd-reports template |
| `R/fcn_plate_layout.R` | `infer_plate_layout()`, copied unchanged from the qmd-reports template |
| `R/fcn_load_data.R` | reads and stacks all plates of a study folder |
| `R/fcn_overview.R` | Overview numbers: progress, positive control, plates table, compound cards and badge rules |
| `R/fcn_setup.R` | Experimental setup: wells of all plates, study design, replicate frames |
| `R/fcn_qc.R` | QC: metrics per run, intensity box statistics, normalization factors, DMSO CVs, PCA, completeness |
| `R/fcn_statistics.R` | Statistics: regulation calls, standard errors, hits per compound, volcano axis limits |
| `R/fcn_proteins.R` | Proteins of interest: hits per protein, abundance rank, minimum detectable change, profiles |
| `R/fcn_analysis.R` | regulation call, fold change to degradation efficiency |
| `R/fcn_plots.R` | plotly plot functions |
| `R/mod_*.R` | one module per tab |
| `R/os_colours.R` | team colours |
| `www/overview.css`, `www/overview.js` | styles and the card wheel |
