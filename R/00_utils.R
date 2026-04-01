# R/00_utils.R
library(dplyr)
library(ggplot2)
library(gt)

options(chromote.timeout = 600)

ensure_chromote_browser <- function(verbose = FALSE) {
  # if already set
  current <- Sys.getenv("CHROMOTE_CHROME", unset = "")
  if (nzchar(current) && file.exists(current)) {
    if (verbose) message("Using CHROMOTE_CHROME: ", current)
    return(invisible(current))
  }
  
  # try to auto-detect via chromote
  chrome <- NULL
  if (requireNamespace("chromote", quietly = TRUE)) {
    chrome <- tryCatch(chromote::find_chrome(), error = function(e) NULL)
  }
  
  if (!is.null(chrome) && nzchar(chrome) && file.exists(chrome)) {
    Sys.setenv(CHROMOTE_CHROME = chrome)
    if (verbose) message("Detected Chrome for chromote: ", chrome)
    return(invisible(chrome))
  }
  
  # fallback common paths
  paths <- c(
    # Windows
    "C:/Program Files/Google/Chrome/Application/chrome.exe",
    "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe",
    "C:/Program Files/Microsoft/Edge/Application/msedge.exe",
    "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe",
    # macOS
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
    # Linux
    "/usr/bin/google-chrome",
    "/usr/bin/google-chrome-stable",
    "/usr/bin/chromium",
    "/usr/bin/chromium-browser",
    "/usr/bin/microsoft-edge"
  )
  
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0) {
    Sys.setenv(CHROMOTE_CHROME = existing[1])
    if (verbose) message("Using fallback browser: ", existing[1])
    return(invisible(existing[1]))
  }
  
  warning(
    "No Chrome/Chromium browser found for chromote. ",
    "gt::gtsave() may fail. Please install Chrome/Edge or set CHROMOTE_CHROME."
  )
  
  invisible(NULL)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_plot_pdf <- function(plot, filepath, width, height, dpi = 300) {
  ensure_dir(dirname(filepath))
  ggsave(filename = filepath, plot = plot, width = width, height = height, dpi = dpi)
  filepath
}

ensure_placeholder_pdf <- function(filepath, message) {
  ensure_dir(dirname(filepath))
  grDevices::pdf(filepath, width = 8, height = 4)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  grid::grid.newpage()
  grid::grid.text(
    message,
    x = 0.5, y = 0.5,
    gp = grid::gpar(fontsize = 12)
  )
  
  invisible(filepath)
}

save_gt_pdf <- function(gt_tbl, filepath, timeout = 60, retries = 2) {
  ensure_dir(dirname(filepath))
  
  old_timeout <- getOption("chromote.timeout")
  on.exit(options(chromote.timeout = old_timeout), add = TRUE)
  options(chromote.timeout = timeout)
  
  last_err <- NULL
  
  for (i in seq_len(retries)) {
    ok <- tryCatch(
      {
        gt::gtsave(gt_tbl, filename = basename(filepath), path = dirname(filepath))
        file.exists(filepath)
      },
      error = function(e) {
        last_err <<- e
        FALSE
      }
    )
    
    if (isTRUE(ok)) {
      return(filepath)
    }
    
    Sys.sleep(1)
  }
  
  warning(
    "gt PDF export failed after retries: ",
    conditionMessage(last_err),
    ". Writing a placeholder PDF instead."
  )
  
  ensure_placeholder_pdf(
    filepath = filepath,
    message = paste0(
      "Table export failed. See console warning for details.\n",
      "This placeholder was written so the pipeline can continue."
    )
  )
}

combine_pdfs <- function(inputs, output) {
  ensure_dir(dirname(output))
  pdftools::pdf_combine(input = inputs, output = output)
  output
}

pretty_virus <- function(x) {
  x <- toupper(x)
  if (x == "EV71") return("EV-A71")
  if (x == "EV68") return("EV-D68")
  x
}
