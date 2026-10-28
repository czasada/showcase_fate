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
    p <- p %>% plotly::add_markers(
      data = present, x = ~plate_col, y = ~plate_row, text = ~text, hoverinfo = "text", name = fill_label,
      marker = c(marker, list(
        color = ~value, colorscale = list(c(0, "#DCE8F5"), c(1, os_colours[1])), showscale = TRUE,
        colorbar = list(title = fill_label, len = 0.8)
      )),
      showlegend = FALSE
    )
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

#' Horizontal bars per sample, coloured by condition.
sample_bar_plot <- function(df, value_col, value_label, sample_order, colours) {
  plotly::plot_ly(
    df,
    x = ~ .data[[value_col]], y = ~sample, color = ~compound, colors = colours, type = "bar", orientation = "h",
    text = ~text, hoverinfo = "text", textposition = "none"
  ) %>%
    finish_plot(xaxis = list(title = value_label), yaxis = category_axis(rev(sample_order)), barmode = "stack")
}

intensity_distribution_plot <- function(intensities, value_col, sample_order, colours) {
  df <- intensities %>% dplyr::filter(!is.na(.data[[value_col]]))
  plotly::plot_ly(
    df,
    x = ~ .data[[value_col]], y = ~sample, color = ~compound, colors = colours, type = "box",
    orientation = "h", boxpoints = FALSE
  ) %>%
    finish_plot(xaxis = list(title = "Log2 intensity"), yaxis = category_axis(rev(sample_order)))
}

cv_plot <- function(cvs, condition_order, colours) {
  medians <- cvs %>% dplyr::summarise(median_cv = stats::median(cv), .by = compound)
  plotly::plot_ly(
    cvs,
    x = ~cv, y = ~compound, color = ~compound, colors = colours, type = "violin", orientation = "h",
    box = list(visible = TRUE), meanline = list(visible = FALSE), points = FALSE, hoverinfo = "x"
  ) %>%
    plotly::add_annotations(
      data = medians, x = 1, xref = "paper", xanchor = "right", y = ~compound,
      text = ~ sprintf("median %.1f %%", median_cv), showarrow = FALSE
    ) %>%
    finish_plot(
      xaxis = list(title = "Coefficient of variation [%]", rangemode = "tozero"),
      yaxis = category_axis(rev(condition_order)), showlegend = FALSE
    )
}

pca_plot <- function(pca, colours, x_pc = "PC1", y_pc = "PC2") {
  i <- as.integer(sub("PC", "", c(x_pc, y_pc)))
  df <- pca$scores %>% dplyr::mutate(text = paste0(sample, "<br>well ", well), replicate = paste0("R", replicate))
  plotly::plot_ly(
    df,
    x = ~ .data[[x_pc]], y = ~ .data[[y_pc]], color = ~compound, colors = colours, symbol = ~replicate,
    type = "scatter", mode = "markers", marker = list(size = 11), text = ~text, hoverinfo = "text"
  ) %>%
    finish_plot(
      xaxis = list(title = sprintf("%s (%.1f %%)", x_pc, pca$variance[i[1]])),
      yaxis = list(title = sprintf("%s (%.1f %%)", y_pc, pca$variance[i[2]]))
    )
}

completeness_plot <- function(intensities, sample_order) {
  df <- intensities %>%
    dplyr::summarise(measured = sum(!is_imputed), imputed = sum(is_imputed), .by = sample) %>%
    tidyr::pivot_longer(c(measured, imputed), names_to = "value_type", values_to = "n") %>%
    dplyr::mutate(text = paste0(sample, "<br>", value_type, ": ", n))
  plotly::plot_ly(
    df,
    x = ~n, y = ~sample, color = ~value_type, colors = c(measured = os_colours[1], imputed = os_colours[2]),
    type = "bar", orientation = "h", text = ~text, hoverinfo = "text", textposition = "none"
  ) %>%
    finish_plot(xaxis = list(title = "Number of proteins"), yaxis = category_axis(rev(sample_order)), barmode = "stack")
}

valid_values_plot <- function(intensities, n_samples) {
  df <- intensities %>%
    dplyr::summarise(n_valid = sum(!is_imputed), .by = protein_id) %>%
    dplyr::count(n_valid) %>%
    dplyr::mutate(text = paste0(n_valid, " of ", n_samples, " samples: ", n, " proteins"))
  plotly::plot_ly(
    df,
    x = ~n_valid, y = ~n, type = "bar", marker = list(color = os_colours[1]),
    text = ~text, hoverinfo = "text", textposition = "none"
  ) %>%
    finish_plot(
      xaxis = list(title = "Samples with a measured value", dtick = 2), yaxis = list(title = "Number of proteins")
    )
}

significant_count_plot <- function(regulated_counts, compound_order) {
  df <- regulated_counts %>%
    dplyr::mutate(n_signed = ifelse(regulation == "down", -n, n), text = paste0(compound, "<br>", regulation, ": ", n))
  plotly::plot_ly(
    df,
    x = ~n_signed, y = ~compound, color = ~regulation, colors = volcano_colours, type = "bar", orientation = "h",
    text = ~n, textposition = "outside", hovertext = ~text, hoverinfo = "text"
  ) %>%
    finish_plot(
      xaxis = list(title = "Number of significant proteins (down | up)", zeroline = TRUE),
      yaxis = category_axis(rev(compound_order)), barmode = "relative"
    )
}

volcano_plot <- function(stats_reg, p_col, p_cutoff, log2fc_cutoff, direction, top_n, highlight_genes, x_range, y_max) {
  p_label <- ifelse(p_col == "p_value_corrected", "BH-adj. p", "p")
  df <- stats_reg %>%
    dplyr::mutate(
      neg_log10_p = -log10(.data[[p_col]]),
      regulation = factor(regulation, levels = c("down", "none", "up")),
      text = paste0(
        "<b>", gene_name, "</b> (", protein_id, ")<br>log2 FC: ", round(log2_fold_change, 2),
        "<br>", p_label, ": ", signif(.data[[p_col]], 2), "<br>imputed values: ", n_imputed
      )
    )
  labelled <- df %>%
    dplyr::filter(regulation != "none") %>%
    dplyr::slice_min(.data[[p_col]], n = top_n, with_ties = FALSE)
  highlight <- df %>% dplyr::filter(gene_name %in% highlight_genes)

  line <- list(color = "#7F7F7F", dash = "dot", width = 1)
  y_cut <- -log10(p_cutoff)
  vline <- function(x) list(type = "line", yref = "paper", y0 = 0, y1 = 1, x0 = x, x1 = x, line = line)
  shapes <- list(list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = y_cut, y1 = y_cut, line = line))
  if (direction != "only_up") shapes <- c(shapes, list(vline(-log2fc_cutoff)))
  if (direction != "only_down") shapes <- c(shapes, list(vline(log2fc_cutoff)))

  p <- plotly::plot_ly() %>%
    plotly::add_markers(
      data = df, x = ~log2_fold_change, y = ~neg_log10_p, color = ~regulation, colors = volcano_colours,
      text = ~text, hoverinfo = "text", type = "scattergl", marker = list(size = 6, opacity = 0.75)
    )
  if (nrow(highlight) > 0) {
    p <- p %>% plotly::add_markers(
      data = highlight, x = ~log2_fold_change, y = ~neg_log10_p, text = ~text, hoverinfo = "text",
      name = "established neosubstrate", type = "scatter", inherit = FALSE,
      marker = list(size = 11, color = "rgba(0,0,0,0)", line = list(color = os_colours[9], width = 2))
    )
  }
  if (nrow(labelled) > 0) {
    p <- p %>% plotly::add_annotations(
      data = labelled, x = ~log2_fold_change, y = ~neg_log10_p, text = ~gene_name,
      showarrow = TRUE, arrowhead = 0, arrowwidth = 0.8, ax = 25, ay = -20, font = list(size = 11)
    )
  }
  finish_plot(
    p,
    shapes = shapes,
    xaxis = list(title = "Log2 fold change", range = x_range, zeroline = FALSE),
    yaxis = list(title = paste0("-log10 (", p_label, ")"), range = c(0, y_max)),
    legend = list(orientation = "h", y = -0.15)
  )
}

overlap_heatmap <- function(overlap, compound_order) {
  m <- overlap %>%
    tidyr::pivot_wider(id_cols = compound_b, names_from = compound_a, values_from = n_shared) %>%
    tibble::column_to_rownames("compound_b") %>%
    as.matrix()
  m <- m[compound_order, compound_order, drop = FALSE]
  plotly::plot_ly(
    x = colnames(m), y = rownames(m), z = m, type = "heatmap",
    colorscale = list(c(0, "#F2F6FB"), c(1, os_colours[2])), showscale = FALSE,
    hovertemplate = "%{x} & %{y}: %{z}<extra></extra>"
  ) %>%
    plotly::add_annotations(
      x = rep(colnames(m), each = nrow(m)), y = rep(rownames(m), times = ncol(m)), text = as.vector(m),
      showarrow = FALSE
    ) %>%
    finish_plot(
      xaxis = category_axis(compound_order, tickangle = -30),
      yaxis = category_axis(rev(compound_order)),
      title = list(text = "Shared hits (same direction)", font = list(size = 14))
    )
}

#' Heatmap of per-sample log2 FC for the selected proteins; rows clustered (complete linkage, Euclidean).
clustermap_plot <- function(fold_changes, protein_ids, proteins, sample_order) {
  m <- fold_changes %>%
    dplyr::filter(protein_id %in% protein_ids) %>%
    tidyr::pivot_wider(id_cols = protein_id, names_from = sample, values_from = log2_fold_change) %>%
    tibble::column_to_rownames("protein_id") %>%
    as.matrix()
  m <- m[, sample_order, drop = FALSE]
  if (nrow(m) > 1) m <- m[stats::hclust(stats::dist(m), method = "complete")$order, , drop = FALSE]
  genes <- proteins$gene_name[match(rownames(m), proteins$protein_id)]
  limit <- max(2, max(abs(m), na.rm = TRUE))
  plotly::plot_ly(
    x = colnames(m), y = genes, z = m, type = "heatmap", zmin = -limit, zmax = limit,
    colorscale = list(c(0, volcano_colours[["down"]]), c(0.5, "#FFFFFF"), c(1, volcano_colours[["up"]])),
    colorbar = list(title = "Log2 FC<br>to mean DMSO"),
    hovertemplate = "%{y}<br>%{x}<br>log2 FC: %{z:.2f}<extra></extra>"
  ) %>%
    finish_plot(
      xaxis = category_axis(sample_order, tickangle = -45),
      yaxis = category_axis(rev(genes), tickfont = list(size = 10))
    )
}

#' Bars (mean of replicates) with replicate points per protein and condition; imputed points open.
#' One small panel per protein, three per row.
poi_bar_plot <- function(fold_changes, value_col, y_label, colours, condition_order, y_range = NULL) {
  genes <- levels(fold_changes$gene_name)
  genes <- genes[genes %in% fold_changes$gene_name]
  panels <- purrr::imap(genes, \(gene, i) {
    df <- fold_changes %>%
      dplyr::filter(gene_name == gene) %>%
      dplyr::mutate(
        value = .data[[value_col]],
        text = paste0(gene, " - ", sample, "<br>", round(value, 2), ifelse(is_imputed, " (imputed)", ""))
      )
    means <- df %>%
      dplyr::summarise(value = mean(value), .by = compound) %>%
      dplyr::mutate(text = paste0(gene, " - ", compound, "<br>mean: ", round(value, 2)))
    plotly::plot_ly(
      means,
      x = ~compound, y = ~value, color = ~compound, colors = colours, type = "bar",
      text = ~text, hoverinfo = "text", textposition = "none", showlegend = FALSE, opacity = 0.85
    ) %>%
      plotly::add_markers(
        data = df, x = ~compound, y = ~value, text = ~text, hoverinfo = "text", showlegend = FALSE, inherit = FALSE,
        marker = list(
          size = 6, color = ifelse(df$is_imputed, "rgba(0,0,0,0)", "#000000"),
          line = list(color = "#000000", width = 1)
        )
      ) %>%
      plotly::layout(
        annotations = list(list(
          text = paste0("<b>", gene, "</b>"), x = 0.5, y = 1.02, xref = "paper", yref = "paper",
          xanchor = "center", yanchor = "bottom", showarrow = FALSE
        )),
        xaxis = category_axis(condition_order, tickangle = -45),
        yaxis = list(title = if ((i - 1) %% 3 == 0) y_label else "", range = y_range)
      )
  })
  n_rows <- ceiling(length(panels) / 3)
  plotly::subplot(
    c(panels, rep(list(plotly::plotly_empty(type = "scatter", mode = "markers")), n_rows * 3 - length(panels))),
    nrows = n_rows, shareX = FALSE, titleY = TRUE, margin = c(0.03, 0.03, 0.06, 0.12)
  ) %>%
    finish_plot()
}
