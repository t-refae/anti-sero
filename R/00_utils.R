# R/00_utils.R
library(dplyr)
library(ggplot2)
library(gt)

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_plot_pdf <- function(plot, filepath, width, height, dpi = 300) {
  ensure_dir(dirname(filepath))
  ggsave(filename = filepath, plot = plot, width = width, height = height, dpi = dpi)
  filepath
}

save_gt_pdf <- function(gt_tbl, filepath) {
  ensure_dir(dirname(filepath))
  gt::gtsave(gt_tbl, filename = basename(filepath), path = dirname(filepath))
  filepath
}

pretty_virus <- function(x) {
  x <- toupper(x)
  if (x == "EV71") return("EV-A71")
  if (x == "EV68") return("EV-D68")
  x
}
