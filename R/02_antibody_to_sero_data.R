# R/02_antibody_to_sero_data.R
library(dplyr)
library(tidyr)
library(stringr)
library(posterior)
library(matrixStats)

# prep antibody raw-well data into Stan list + per-serum metadata 
prep_antibody_data <- function(virus, raw_prefix, data_dir = "data") {
  csv_path <- file.path(data_dir, paste0(raw_prefix, "_raw_well_observations.csv"))
  dat <- read.csv(csv_path) %>%
    mutate(serumID = as.character(serumID)) %>%
    arrange(serumID, dilutions)
  
  # replicate columns (rep1, rep2, ...)
  n_reps <- length(grep("^rep", names(dat)))
  stopifnot(n_reps > 0)
  stopifnot(max(dat$outcome) <= n_reps)
  
  d_levels <- sort(unique(dat$dilutions))
  
  z_mat <- dat %>%
    mutate(dilutions = factor(dilutions, levels = d_levels)) %>%
    select(serumID, dilutions, outcome) %>%
    pivot_wider(names_from = dilutions, values_from = outcome) %>%
    arrange(serumID) %>%
    select(-serumID) %>%
    as.matrix()
  
  serum_meta <- dat %>%
    group_by(serumID) %>%
    summarise(
      age  = suppressWarnings(as.numeric(first(Age))),
      year = suppressWarnings(as.numeric(first(Year))),
      .groups = "drop"
    ) %>%
    arrange(serumID)
  
  stopifnot(nrow(serum_meta) == nrow(z_mat))
  
  stan_data <- list(
    n_individuals = nrow(z_mat),
    n_dilutions   = length(d_levels),
    n_replicates  = n_reps,
    d             = d_levels,
    z             = z_mat
  )
  
  list(
    virus = virus,
    raw_prefix = raw_prefix,
    csv_path = csv_path,
    stan_data = stan_data,
    serum_meta = serum_meta
  )
}

# extract log(phi) draws from antibody fit (optionally thin for speed)
extract_log_phi_draws <- function(draws_obj, n_individuals, thin_to = 500, seed = 1) {
  ddf <- posterior::as_draws_df(draws_obj)
  
  phi_cols <- grep("^phi\\[", names(ddf), value = TRUE)
  idx <- as.integer(stringr::str_match(phi_cols, "^phi\\[(\\d+)\\]")[, 2])
  phi_cols <- phi_cols[order(idx)]
  
  stopifnot(length(phi_cols) == n_individuals)
  
  phi_mat <- as.matrix(ddf[, phi_cols, drop = FALSE]) # [draw x person]
  if (!is.null(thin_to) && thin_to < nrow(phi_mat)) {
    set.seed(seed)
    keep <- sample.int(nrow(phi_mat), thin_to)
    phi_mat <- phi_mat[keep, , drop = FALSE]
  }
  
  log(phi_mat)
}

# stan data for serocatalytic model
make_sero_stan_data <- function(serum_meta, log_phi_draws) {
  stopifnot(ncol(log_phi_draws) == nrow(serum_meta))
  
  ages <- as.integer(round(serum_meta$age))
  age_max <- max(ages, na.rm = TRUE)
  
  list(
    n = nrow(serum_meta),
    ages = ages,
    age_max = as.integer(age_max),
    n_phi_draws = as.integer(nrow(log_phi_draws)),
    log_phi_draws = log_phi_draws
  )
}

# plotting data (phi medians + ages)
make_sero_plot_data <- function(serum_meta, log_phi_draws) {
  ages <- as.integer(round(serum_meta$age))
  phi_median <- exp(matrixStats::colMedians(log_phi_draws))
  
  list(
    n = length(ages),
    ages = ages,
    age_max = max(ages, na.rm = TRUE),
    phi_median = phi_median
  )
}