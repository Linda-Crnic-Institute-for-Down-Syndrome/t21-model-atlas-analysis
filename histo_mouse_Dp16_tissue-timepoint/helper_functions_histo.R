# Histopathology analysis plotting and export helpers
#
# Plot theme ----
ggplot2::theme_set(
  ggplot2::theme_gray(base_size = 12, base_family = "Arial") +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = "transparent"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      axis.text = ggplot2::element_text(color = "black", size = 14),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank()
    )
)
#

# Excel export ----
export_excel <- function(named_list, filename = "") {
  wb <- openxlsx::createWorkbook()
  Map(function(data, name) {
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeData(wb, name, data)
  }, named_list, names(named_list))
  path <- here("results", paste0(out_file_prefix, filename, ".xlsx"))
  openxlsx::saveWorkbook(wb, file = path, overwrite = TRUE)
  cat("Saved as:", path, "\n")
}
#

# Individual score plot ----
# Jittered points per animal with a median bar, one facet per panel label.
# x: tibble with Genotype, Score_value and the facet column(s)
individual_score_plot <- function(x, facets, title = "", y_label = "Score", ncol = NULL, colors = standard_colors) {
  x |>
    ggplot(aes(x = Genotype, y = Score_value)) +
    geom_jitter(width = 0.2, height = 0, aes(color = Genotype)) +
    stat_boxplot(
      aes(ymin = after_stat(middle), ymax = after_stat(middle)),
      geom = "errorbar", linewidth = 0.4, width = 0.6, colour = "black"
    ) +
    scale_color_manual(values = colors) +
    scale_y_continuous(breaks = \(x) seq(0, max(x), by = 1),
                       expand = expansion(mult = c(0.05, 0.15))) +
    facet_wrap(facets, scales = "free_y", ncol = ncol) +
    theme(
      legend.position = "none",
      aspect.ratio = 1.3,
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold", size = 10),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = title, x = "Genotype", y = y_label)
}
#
