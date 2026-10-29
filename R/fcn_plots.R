# Plot functions for the dashboard, built with native plotly (no ggplotly, so no dependency on
# ggplot2/plotly version pairs).

finish_plot <- function(p, ...) {
  p %>%
    plotly::layout(font = list(size = 13), ...) %>%
    plotly::config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d"))
}

category_axis <- function(levels, title = "", ...) {
  list(title = title, type = "category", categoryorder = "array", categoryarray = levels, ...)
}

#' Plate map: one square per well, coloured by a discrete column (named colours) or a continuous one.
plate_plot <- function(wells, fill_col, fill_label, n_rows, n_cols, colours = NULL, marker_size = 20) {
  rows <- rev(LETTERS[seq_len(n_rows)])
  wells <- wells %>% dplyr::mutate(value = .data[[fill_col]])
  marker <- list(symbol = "square", size = marker_size, line = list(color = "#FFFFFF", width = 1))
  p <- plotly::plot_ly()
  if (is.null(colours)) {
    missing <- wells %>% dplyr::filter(is.na(value))
    present <- wells %>% dplyr::filter(!is.na(value))
    if (nrow(missing) > 0) {
      p <- p %>% plotly::add_markers(
        data = missing, x = ~plate_col, y = ~plate_row, text = ~text, hoverinfo = "text", name = "n/a",
        marker = c(marker, list(color = "#F2F2F2")), showlegend = FALSE
      )
    }
    if (nrow(present) > 0) {
      p <- p %>% plotly::add_markers(
        data = present, x = ~plate_col, y = ~plate_row, text = ~text, hoverinfo = "text", name = fill_label,
        marker = c(marker, list(
          color = ~value, colorscale = list(c(0, "#DCE8F5"), c(1, os_colours[1])), showscale = TRUE,
          colorbar = list(title = fill_label, len = 0.8)
        )),
        showlegend = FALSE
      )
    }
  } else {
    for (level in names(colours)) {
      level_wells <- wells %>% dplyr::filter(as.character(value) == level)
      if (nrow(level_wells) == 0) next
      p <- p %>% plotly::add_markers(
        data = level_wells, x = ~plate_col, y = ~plate_row, text = ~text, hoverinfo = "text", name = level,
        marker = c(marker, list(color = colours[[level]]))
      )
    }
    empty <- wells %>% dplyr::filter(is.na(value))
    if (nrow(empty) > 0) {
      p <- p %>% plotly::add_markers(
        data = empty, x = ~plate_col, y = ~plate_row, text = ~text, hoverinfo = "text", name = "empty",
        marker = c(marker, list(color = "#F2F2F2"))
      )
    }
  }
  finish_plot(
    p,
    xaxis = list(title = "", side = "top", dtick = 1, range = c(0.5, n_cols + 0.5), showgrid = FALSE, zeroline = FALSE),
    yaxis = category_axis(rows, showgrid = FALSE),
    legend = list(title = list(text = fill_label))
  )
}

#' Plate blocks for the QC strips: shaded background of every second plate and a label per plate.
plate_block_layout <- function(runs) {
  plates <- runs %>%
    dplyr::summarise(x0 = min(x) - 0.5, x1 = max(x) + 0.5, .by = plate_id) %>%
    dplyr::mutate(i = dplyr::row_number())
  list(
    shapes = purrr::pmap(plates %>% dplyr::filter(i %% 2 == 0), \(plate_id, x0, x1, i) {
      list(type = "rect", x0 = x0, x1 = x1, y0 = 0, y1 = 1, yref = "paper", fillcolor = "#f4f6f8",
           line = list(width = 0), layer = "below")
    }),
    annotations = purrr::pmap(plates, \(plate_id, x0, x1, i) {
      list(x = (x0 + x1) / 2, y = 1, yref = "paper", yanchor = "bottom", showarrow = FALSE,
           text = paste("Plate", sub("P", "", plate_id)))
    })
  )
}

#' QC strip: one dot per run in measurement order over all plates, flagged runs highlighted, optional threshold
#' line per plate (share of the plate median).
run_metric_plot <- function(runs, metric, label, threshold_share = NULL) {
  runs <- runs %>%
    dplyr::mutate(
      value = .data[[metric]],
      colour = ifelse(qc_status == "flagged", "#C77E00", well_type_colours[well_type]),
      text = sprintf(
        "<b>%s</b> (plate %s, %s)<br>%s: %s<br>%.0f %% of plate median<br>QC: %s", sample, sub("P", "", plate_id),
        position, label, format(signif(value, 4), big.mark = ","), value / stats::median(value) * 100, qc_status
      ),
      .by = plate_id
    )
  blocks <- plate_block_layout(runs)
  shapes <- blocks$shapes
  if (!is.null(threshold_share)) {
    thresholds <- runs %>%
      dplyr::summarise(
        y = threshold_share * stats::median(value), x0 = min(x) - 0.5, x1 = max(x) + 0.5,
        .by = plate_id
      )
    shapes <- c(shapes, purrr::pmap(thresholds, \(plate_id, y, x0, x1) {
      list(type = "line", x0 = x0, x1 = x1, y0 = y, y1 = y, line = list(color = "#C77E00", dash = "dash", width = 1))
    }))
  }
  plotly::plot_ly(
    runs,
    x = ~x, y = ~value, type = "scatter", mode = "markers", text = ~text, hoverinfo = "text",
    marker = list(color = ~colour, size = 7)
  ) %>%
    finish_plot(
      xaxis = list(title = "Run (measurement order)", showticklabels = FALSE, zeroline = FALSE),
      yaxis = list(title = label), shapes = shapes, annotations = blocks$annotations, showlegend = FALSE
    )
}

#' Intensity distribution per run from precomputed box statistics.
intensity_box_plot <- function(box_stats, label) {
  plotly::plot_ly(
    type = "box", x = box_stats$sample, q1 = box_stats$q1, median = box_stats$median, q3 = box_stats$q3,
    lowerfence = box_stats$lower, upperfence = box_stats$upper,
    fillcolor = os_colours[2], line = list(color = os_colours[1], width = 1), hoverinfo = "y"
  ) %>%
    finish_plot(
      xaxis = category_axis(box_stats$sample, tickangle = -90, tickfont = list(size = 8)),
      yaxis = list(title = label), showlegend = FALSE
    )
}

#' Normalization factor per run in measurement order over all plates.
normalization_factor_plot <- function(factors) {
  blocks <- plate_block_layout(factors)
  factors <- factors %>%
    dplyr::mutate(text = sprintf("<b>%s</b> (plate %s)<br>factor: %.2f", sample, sub("P", "", plate_id), factor))
  plotly::plot_ly(
    factors,
    x = ~x, y = ~factor, type = "scatter", mode = "markers", text = ~text, hoverinfo = "text",
    marker = list(color = unname(well_type_colours[factors$well_type]), size = 7)
  ) %>%
    finish_plot(
      xaxis = list(title = "Run (measurement order)", showticklabels = FALSE, zeroline = FALSE),
      yaxis = list(title = "Normalization factor (log2)"), shapes = blocks$shapes,
      annotations = blocks$annotations, showlegend = FALSE
    )
}

#' CV distributions of the DMSO wells per plate, raw and normalized, with a 20 % reference line.
dmso_cv_plot <- function(cvs) {
  cvs <- cvs %>%
    dplyr::mutate(
      type = factor(type, levels = c("raw", "normalized")),
      box = paste0("Plate ", sub("P", "", plate_id), " · ", type)
    ) %>%
    dplyr::arrange(plate_id, type)
  plotly::plot_ly(
    cvs,
    x = ~box, y = ~cv, color = ~type, colors = c(raw = os_colours[2], normalized = os_colours[1]),
    type = "box", boxpoints = FALSE
  ) %>%
    finish_plot(
      xaxis = category_axis(unique(cvs$box)), yaxis = list(title = "CV (%)", rangemode = "tozero"),
      shapes = list(list(
        type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 20, y1 = 20,
        line = list(color = "#a6a6a6", dash = "dash", width = 1)
      ))
    )
}

#' PCA scores of all runs, coloured by plate or well type.
qc_pca_plot <- function(pca, colour_by, x_pc = "PC1", y_pc = "PC2") {
  scores <- pca$scores %>%
    dplyr::mutate(
      group = if (colour_by == "plate_id") paste("Plate", sub("P", "", plate_id)) else well_type,
      text = sprintf("<b>%s</b> (plate %s)<br>QC: %s", sample, sub("P", "", plate_id), qc_status)
    )
  colours <- if (colour_by == "plate_id") {
    stats::setNames(os_colours_selected[seq_along(unique(scores$group))], sort(unique(scores$group)))
  } else {
    well_type_colours
  }
  axis_title <- \(pc) sprintf("%s (%.1f %%)", pc, pca$variance[as.integer(sub("PC", "", pc))])
  plotly::plot_ly(
    scores,
    x = scores[[x_pc]], y = scores[[y_pc]], color = ~group, colors = colours, text = ~text, hoverinfo = "text",
    type = "scatter", mode = "markers", marker = list(size = 8)
  ) %>%
    finish_plot(xaxis = list(title = axis_title(x_pc)), yaxis = list(title = axis_title(y_pc)))
}

scree_plot <- function(variance, n = 10) {
  k <- seq_len(min(n, length(variance)))
  df <- tibble::tibble(pc = paste0("PC", k), variance = variance[k])
  plotly::plot_ly(df, x = ~pc, y = ~variance, type = "bar", marker = list(color = os_colours[1]),
                  hovertemplate = "%{x}: %{y:.1f} %<extra></extra>") %>%
    finish_plot(xaxis = category_axis(df$pc), yaxis = list(title = "Explained variance (%)"))
}

#' Number of proteins by how many runs of a plate measured them.
valid_values_plot <- function(per_protein) {
  per_protein <- per_protein %>% dplyr::mutate(plate = paste("Plate", sub("P", "", plate_id)))
  plotly::plot_ly(
    per_protein,
    x = ~valid, y = ~n, color = ~plate, colors = os_colours_selected[seq_along(unique(per_protein$plate))],
    type = "bar", hovertemplate = "%{y} proteins measured in %{x} runs<extra></extra>"
  ) %>%
    finish_plot(
      barmode = "group", xaxis = list(title = "Runs with a value"), yaxis = list(title = "Proteins", type = "log")
    )
}

#' Statistics summary: hits per compound, sorted, down hits split into known neosubstrates and other proteins.
#' Clicking a bar reports the compound through event_data("plotly_click", source).
hit_count_plot <- function(counts, source) {
  long <- counts %>%
    dplyr::transmute(compound, plate_id,
      `known neosubstrates` = neo_down, `other proteins down` = down - neo_down, `up` = up
    ) %>%
    tidyr::pivot_longer(-c(compound, plate_id), names_to = "kind", values_to = "n") %>%
    dplyr::filter(n > 0) %>%
    dplyr::mutate(text = sprintf("<b>%s</b> (plate %s)<br>%s: %d", compound, sub("P", "", plate_id), kind, n))
  plotly::plot_ly(
    long,
    x = ~compound, y = ~n, color = ~kind, type = "bar", text = ~text, hoverinfo = "text", textposition = "none",
    customdata = ~compound, source = source,
    colors = c(`known neosubstrates` = "#c40a4f", `other proteins down` = "#C8C8C8", up = os_colours[1])
  ) %>%
    finish_plot(
      barmode = "stack", xaxis = category_axis(counts$compound, tickangle = -90, tickfont = list(size = 9)),
      yaxis = list(title = "Proteins"), legend = list(orientation = "h", x = 0, y = 1.12)
    ) %>%
    plotly::event_register("plotly_click")
}

#' Volcano of one compound with the -log10 p-value (y_mode "p") or the -log10 standard error (y_mode "se") on the y
#' axis, so precise measurements sit high in both (SE volcano after the 0609 volcano report). Points beyond the x
#' limit are drawn on the edge. Clicking a point reports its scout_id through event_data("plotly_click", source).
volcano_plot <- function(reg, y_mode, p_col, p_cutoff, log2fc_cutoff, direction, top_n, highlight_genes,
                         source = NULL, selected_id = NULL) {
  x_limit <- volcano_x_limit(reg$log2_fold_change)
  p_label <- ifelse(p_col == "p_value_corrected", "BH-adj. p", "p")
  df <- reg %>%
    dplyr::mutate(
      x = pmax(pmin(log2_fold_change, x_limit), -x_limit),
      y = if (y_mode == "p") -log10(.data[[p_col]]) else -log10(standard_error),
      regulation = factor(regulation, levels = c("down", "none", "up")),
      text = sprintf(
        "<b>%s</b><br>log2 FC: %.2f<br>%s: %s<br>SE: %.3f<br>valid values: %s", gene_name, log2_fold_change,
        p_label, signif(.data[[p_col]], 2), standard_error, valid_value_string
      )
    ) %>%
    dplyr::filter(!is.na(y))
  y_range <- if (y_mode == "p") {
    c(0, -log10(min(min(reg[[p_col]], na.rm = TRUE), volcano_limits$min_p) * 0.1))
  } else {
    -log10(c(
      max(max(reg$standard_error, na.rm = TRUE), volcano_limits$max_standard_error),
      min(min(reg$standard_error, na.rm = TRUE), volcano_limits$min_standard_error)
    ))
  }
  labelled <- df %>%
    dplyr::filter(regulation != "none") %>%
    dplyr::slice_min(.data[[p_col]], n = top_n, with_ties = FALSE)
  highlight <- df %>% dplyr::filter(gene_name %in% highlight_genes)

  line <- list(color = "#7F7F7F", dash = "dot", width = 1)
  vline <- function(x) list(type = "line", yref = "paper", y0 = 0, y1 = 1, x0 = x, x1 = x, line = line)
  shapes <- list()
  if (direction != "only_up") shapes <- c(shapes, list(vline(-log2fc_cutoff)))
  if (direction != "only_down") shapes <- c(shapes, list(vline(log2fc_cutoff)))
  if (y_mode == "p") {
    y_cut <- -log10(p_cutoff)
    shapes <- c(shapes, list(list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = y_cut, y1 = y_cut, line = line)))
  }

  p <- plotly::plot_ly(source = source) %>%
    plotly::add_markers(
      data = df, x = ~x, y = ~y, color = ~regulation, colors = volcano_colours, customdata = ~scout_id,
      text = ~text, hoverinfo = "text", type = "scattergl", marker = list(size = 6, opacity = 0.75)
    )
  selected <- df %>% dplyr::filter(scout_id %in% selected_id)
  if (nrow(selected) > 0) {
    p <- p %>% plotly::add_markers(
      data = selected, x = ~x, y = ~y, text = ~text, hoverinfo = "text", customdata = ~scout_id,
      name = "selected", type = "scatter", inherit = FALSE,
      marker = list(size = 15, color = "rgba(0,0,0,0)", line = list(color = "#000000", width = 2))
    )
  }
  if (nrow(highlight) > 0) {
    p <- p %>% plotly::add_markers(
      data = highlight, x = ~x, y = ~y, text = ~text, hoverinfo = "text", customdata = ~scout_id,
      name = "known neosubstrate", type = "scatter", inherit = FALSE,
      marker = list(size = 11, color = "rgba(0,0,0,0)", line = list(color = os_colours[9], width = 2))
    )
  }
  if (nrow(labelled) > 0) {
    p <- p %>% plotly::add_annotations(
      data = labelled, x = ~x, y = ~y, text = ~gene_name,
      showarrow = TRUE, arrowhead = 0, arrowwidth = 0.8, ax = 25, ay = -20, font = list(size = 11)
    )
  }
  finish_plot(
    p,
    shapes = shapes,
    xaxis = list(title = "Log2 fold change", range = c(-x_limit, x_limit) * 1.1, zeroline = FALSE),
    yaxis = list(
      title = if (y_mode == "p") paste0("-log10 (", p_label, ")") else "-log10 (standard error)", range = y_range
    ),
    legend = list(orientation = "h", y = -0.18)
  ) %>%
    plotly::event_register("plotly_click")
}

#' Replicate values of one protein: DMSO wells of the plate next to the compound's wells, per compound shown.
#' Imputed values are drawn hollow. The first compound is the volcano's compound, further ones (magic toggle) follow.
protein_replicate_plot <- function(values, gene) {
  values <- values %>%
    dplyr::mutate(
      colour = dplyr::case_when(
        role == "DMSO" ~ vehicle_colour, role == "compound" ~ "#c40a4f", .default = os_colours[1]
      ),
      symbol = ifelse(is_imputed, "circle-open", "circle"),
      text = sprintf(
        "<b>%s</b> %s (plate %s)<br>log2 intensity: %.2f%s", group, replicate, sub("P", "", plate_id), value,
        ifelse(is_imputed, "<br>imputed", "")
      )
    )
  means <- values %>% dplyr::summarise(mean = mean(value), .by = group)
  groups <- unique(values$group)
  plotly::plot_ly(
    values,
    x = ~group, y = ~value, type = "scatter", mode = "markers", text = ~text, hoverinfo = "text",
    marker = list(color = ~colour, symbol = ~symbol, size = 11, line = list(color = ~colour, width = 2))
  ) %>%
    finish_plot(
      title = list(text = gene, x = 0.02),
      xaxis = category_axis(groups, tickangle = -30), yaxis = list(title = "Normalized log2 intensity"),
      showlegend = FALSE,
      shapes = purrr::map2(match(means$group, groups) - 1, means$mean, \(i, m) {
        list(type = "line", x0 = i - 0.3, x1 = i + 0.3, y0 = m, y1 = m, line = list(color = "#444444", width = 2))
      })
    )
}

# Hover text of one compound in the protein views.
protein_profile_text <- function(profile) {
  sprintf(
    "<b>%s</b> (plate %s)<br>log2 FC: %.2f (%.0f %% reduction)<br>BH-adj. p: %s<br>degrades %d proteins in total",
    profile$compound, sub("P", "", profile$plate_id), profile$log2_fold_change,
    fold_change_to_degradation_efficiency(2^profile$log2_fold_change), signif(profile$p_value_corrected, 2),
    profile$n_degraded
  )
}

# Reference lines of the protein views on the chosen metric: hit cutoff and minimum detectable change.
protein_reference_lines <- function(metric, mdc, axis = "y") {
  lines <- list(
    list(value = protein_metric_value(hit_log2fc_cutoff, metric), label = "hit cutoff (50 %)", colour = "#c40a4f"),
    list(value = protein_metric_value(mdc, metric), label = "min. detectable change", colour = "#7F7F7F")
  )
  other <- if (axis == "y") "x" else "y"
  list(
    shapes = lapply(lines, \(l) {
      shape <- list(type = "line", line = list(color = l$colour, dash = "dash", width = 1))
      shape[[paste0(other, "ref")]] <- "paper"
      shape[[paste0(other, "0")]] <- 0
      shape[[paste0(other, "1")]] <- 1
      shape[[paste0(axis, "0")]] <- l$value
      shape[[paste0(axis, "1")]] <- l$value
      shape
    }),
    annotations = lapply(lines, \(l) {
      a <- list(text = l$label, showarrow = FALSE, font = list(size = 10, color = l$colour))
      a[[paste0(other, "ref")]] <- "paper"
      a[[other]] <- 1
      a[[axis]] <- l$value
      anchors <- if (axis == "y") {
        list(xanchor = "right", yanchor = "bottom")
      } else {
        list(yanchor = "top", xanchor = "left")
      }
      c(a, anchors)
    })
  )
}

#' Proteins of interest, bubble view: one bubble per compound in ID order (plates shaded), y = metric, filled = hit,
#' size = proteins the compound degrades in total (small filled bubbles are selective degraders).
protein_bubble_plot <- function(profile, metric, mdc) {
  profile <- profile %>%
    dplyr::mutate(
      value = protein_metric_value(log2_fold_change, metric), text = protein_profile_text(profile),
      colour = ifelse(hit, "rgba(196,10,79,0.8)", "rgba(166,166,166,0.35)")
    )
  blocks <- plate_block_layout(profile %>% dplyr::mutate(x = order))
  refs <- protein_reference_lines(metric, mdc)
  plotly::plot_ly(
    profile,
    x = ~order, y = ~value, type = "scatter", mode = "markers", text = ~text, hoverinfo = "text",
    marker = list(
      color = ~colour, size = ~n_degraded, sizemode = "area", sizeref = max(profile$n_degraded, 1) / 30^2,
      sizemin = 4, line = list(width = 0)
    )
  ) %>%
    finish_plot(
      xaxis = list(title = sprintf("Compounds in ID order (%d with results)", nrow(profile)), showticklabels = FALSE,
                   zeroline = FALSE),
      yaxis = list(title = protein_metric_label[[metric]]),
      shapes = c(blocks$shapes, refs$shapes), annotations = c(blocks$annotations, refs$annotations),
      showlegend = FALSE
    )
}

#' Proteins of interest, volcano view: this protein in every compound, log2 FC vs -log10 BH-adjusted p.
protein_volcano_plot <- function(profile, top_n = 5) {
  profile <- profile %>%
    dplyr::mutate(
      y = -log10(p_value_corrected), text = protein_profile_text(profile),
      regulation = factor(ifelse(hit, "down", "none"), levels = c("down", "none"))
    )
  labelled <- profile %>%
    dplyr::filter(hit) %>%
    dplyr::slice_min(log2_fold_change, n = top_n, with_ties = FALSE)
  line <- list(color = "#7F7F7F", dash = "dot", width = 1)
  p <- plotly::plot_ly(
    profile,
    x = ~log2_fold_change, y = ~y, color = ~regulation, colors = volcano_colours, text = ~text, hoverinfo = "text",
    type = "scatter", mode = "markers", marker = list(size = 9, opacity = 0.8)
  )
  if (nrow(labelled) > 0) {
    p <- p %>% plotly::add_annotations(
      data = labelled, x = ~log2_fold_change, y = ~y, text = ~compound,
      showarrow = TRUE, arrowhead = 0, arrowwidth = 0.8, ax = 30, ay = -20, font = list(size = 11)
    )
  }
  finish_plot(
    p,
    xaxis = list(title = "Log2 fold change", zeroline = FALSE),
    yaxis = list(title = "-log10 (BH-adj. p)", rangemode = "tozero"),
    shapes = list(
      list(type = "line", yref = "paper", y0 = 0, y1 = 1, x0 = hit_log2fc_cutoff, x1 = hit_log2fc_cutoff, line = line),
      list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = -log10(hit_p_cutoff), y1 = -log10(hit_p_cutoff),
           line = line)
    ),
    legend = list(orientation = "h", y = -0.18)
  )
}

#' Proteins of interest, ranked view: compounds by reduction of this protein (top n), hits coloured.
protein_ranked_plot <- function(profile, metric, mdc, n = 25) {
  ranked <- profile %>% dplyr::slice_min(log2_fold_change, n = n, with_ties = FALSE)
  ranked <- ranked %>%
    dplyr::mutate(value = protein_metric_value(log2_fold_change, metric), text = protein_profile_text(ranked))
  refs <- protein_reference_lines(metric, mdc, axis = "x")
  plotly::plot_ly(
    ranked,
    x = ~value, y = ~compound, type = "bar", orientation = "h", text = ~text, hoverinfo = "text",
    textposition = "none", marker = list(color = ifelse(ranked$hit, "#c40a4f", "#C8C8C8"))
  ) %>%
    finish_plot(
      xaxis = list(title = protein_metric_label[[metric]]),
      yaxis = category_axis(rev(ranked$compound), tickfont = list(size = 10)),
      shapes = refs$shapes, annotations = refs$annotations, showlegend = FALSE
    )
}

#' Overview: cumulative runs acquired against plan, with a marker for today.
runs_plan_plot <- function(series, today = Sys.Date()) {
  acquired <- series$acquired
  last <- acquired[nrow(acquired), ]
  plotly::plot_ly() %>%
    plotly::add_lines(
      data = series$plan, x = ~date, y = ~runs, name = "plan",
      line = list(shape = "hv", dash = "dash", color = os_colours[2]),
      hovertemplate = "%{x|%d %b}: %{y} runs planned<extra></extra>"
    ) %>%
    plotly::add_lines(
      data = acquired, x = ~date, y = ~runs, name = "acquired",
      line = list(shape = "hv", width = 3, color = os_colours[1]),
      hovertemplate = "%{x|%d %b}: %{y} runs acquired<extra></extra>"
    ) %>%
    plotly::add_markers(
      data = last, x = ~date, y = ~runs, marker = list(color = os_colours[1], size = 9), showlegend = FALSE,
      hoverinfo = "skip"
    ) %>%
    finish_plot(
      xaxis = list(title = "", tickformat = "%d %b"), yaxis = list(title = "Runs", rangemode = "tozero"),
      legend = list(orientation = "h", x = 0, y = 1.12),
      shapes = list(list(
        type = "line", x0 = today, x1 = today, y0 = 0, y1 = 1, yref = "paper",
        line = list(color = os_colours[3], dash = "dot", width = 1)
      )),
      annotations = list(
        list(x = today, y = 1, yref = "paper", text = "today", showarrow = FALSE, xanchor = "left",
             font = list(color = os_colours[3], size = 11)),
        list(x = last$date, y = last$runs, text = paste(last$runs, "runs"), showarrow = FALSE, yshift = 14)
      )
    )
}

#' Overview: positive control reduction per plate and replicate, acceptance band, hollow markers for plates
#' without results.
positive_control_plot <- function(replicates, plates, acceptance = pos_ctrl_acceptance) {
  missing <- setdiff(plates, replicates$plate_id)
  replicates <- replicates %>%
    dplyr::mutate(text = sprintf("Plate %s, %s: %.0f %% reduction", sub("P", "", plate_id), replicate, reduction))
  p <- plotly::plot_ly() %>%
    plotly::add_markers(
      data = replicates, x = ~plate_id, y = ~reduction, text = ~text, hoverinfo = "text",
      marker = list(color = ifelse(replicates$reduction >= acceptance, os_colours[1], os_colours[3]), size = 12),
      showlegend = FALSE
    )
  if (length(missing) > 0) {
    p <- p %>% plotly::add_markers(
      x = missing, y = rep(acceptance - 15, length(missing)), text = "not yet acquired", hoverinfo = "text",
      marker = list(symbol = "circle-open", color = "#a6a6a6", size = 12), showlegend = FALSE
    )
  }
  finish_plot(
    p,
    xaxis = category_axis(plates, ticktext = paste("Plate", sub("P", "", plates)), tickvals = plates),
    yaxis = list(title = "GSPT1 reduction (%)", range = c(50, 100), ticksuffix = " %"),
    shapes = list(
      list(type = "rect", x0 = 0, x1 = 1, xref = "paper", y0 = acceptance, y1 = 100,
           fillcolor = "#eaf3de", opacity = 0.6, line = list(width = 0), layer = "below"),
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = acceptance, y1 = acceptance,
           line = list(color = "#6C7A05", dash = "dash", width = 1))
    ),
    annotations = list(list(
      x = 1, xref = "paper", y = acceptance, text = paste0("acceptance ", acceptance, " %"), showarrow = FALSE,
      xanchor = "right", yanchor = "top", font = list(color = "#6C7A05", size = 11)
    ))
  )
}
