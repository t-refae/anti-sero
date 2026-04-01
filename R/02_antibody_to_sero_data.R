# R/02_antibody_to_sero_data.R
library(dplyr)
library(tidyr)
library(stringr)
library(posterior)
library(matrixStats)

# prep antibody raw-well data into Stan list + per-serum metadata 
prep_antibody_data <- function(virus, raw_prefix, data_dir = "data/raw") {
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

#### Reed-Muench ####
make_reed_muench_data <- function() {
  raw_df <- data.frame(
    serumID = rep("RM1938", 9),
    dilutions = c(1, 2, 4, 8, 16, 32, 64, 128, 256),
    outcome = c(6, 6, 5, 6, 4, 2, 2, 0, 1),
    n_replicates = rep(6L, 9),
    stringsAsFactors = FALSE
  )
  
  stan_data <- list(
    n_individuals = 1L,
    n_dilutions = nrow(raw_df),
    n_replicates = 6L,
    d = raw_df$dilutions,
    z = matrix(raw_df$outcome, nrow = 1)
  )
  
  list(raw_df = raw_df, stan_data = stan_data)
}

reed_muench_endpoint <- function(raw_df) {
  df <- raw_df %>%
    dplyr::arrange(dilutions) %>%
    dplyr::mutate(
      dead = n_replicates - outcome,
      cum_survival = rev(cumsum(rev(outcome))),
      cum_deaths = cumsum(dead),
      pct_mortality = 100 * cum_deaths / (cum_deaths + cum_survival)
    )
  
  below <- max(which(df$pct_mortality < 50))
  above <- min(which(df$pct_mortality > 50))
  
  prop_dist <- (50 - df$pct_mortality[below]) /
    (df$pct_mortality[above] - df$pct_mortality[below])
  
  10^(
    log10(df$dilutions[below]) +
      prop_dist * log10(df$dilutions[above] / df$dilutions[below])
  )
}

summarise_reed_muench_fit <- function(draws_obj) {
  ddf <- posterior::as_draws_df(draws_obj)
  phi <- ddf[["phi[1]"]]
  k1  <- ddf[["k1"]]
  
  out <- rbind(
    data.frame(
      parameter = "phi",
      median = stats::median(phi),
      q2.5 = unname(stats::quantile(phi, 0.025)),
      q97.5 = unname(stats::quantile(phi, 0.975))
    ),
    data.frame(
      parameter = "k1",
      median = stats::median(k1),
      q2.5 = unname(stats::quantile(k1, 0.025)),
      q97.5 = unname(stats::quantile(k1, 0.975))
    ),
    data.frame(
      parameter = "endpoint_50",
      median = stats::median(phi),
      q2.5 = unname(stats::quantile(phi, 0.025)),
      q97.5 = unname(stats::quantile(phi, 0.975))
    )
  )
  
  rownames(out) <- NULL
  out
}

#### LOOCV ####
extract_log_likelihood_draws <- function(draws_obj, var = "log_likelihood") {
  ddf <- posterior::as_draws_df(draws_obj)
  
  ll_cols <- grep(paste0("^", var, "\\["), names(ddf), value = TRUE)
  stopifnot(length(ll_cols) > 0)
  
  idx <- as.integer(stringr::str_match(
    ll_cols,
    paste0("^", var, "\\[(\\d+)\\]$")
  )[, 2])
  
  ll_cols <- ll_cols[order(idx)]
  as.matrix(ddf[, ll_cols, drop = FALSE])
}

compare_si_sii_loo <- function(si_draws, sii_draws, virus) {
  loo_si  <- loo::loo(extract_log_likelihood_draws(si_draws))
  loo_sii <- loo::loo(extract_log_likelihood_draws(sii_draws))
  
  cmp <- loo::loo_compare(list(SI = loo_si, SII = loo_sii))
  
  out <- data.frame(
    Virus = virus,
    Model = rownames(cmp),
    ELPD_diff = unname(cmp[, "elpd_diff"]),
    SE_diff = unname(cmp[, "se_diff"]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  out
}

combine_loo_tables <- function(ev71_loo, cva6_loo, ev68_loo) {
  dplyr::bind_rows(cva6_loo, ev71_loo, ev68_loo) %>%
    dplyr::mutate(
      Virus = factor(Virus, levels = c("CVA6", "EV71", "EV68")),
      Model = factor(Model, levels = c("SII", "SI"))
    ) %>%
    dplyr::arrange(Virus, dplyr::desc(Model == "SII"), Model) %>%
    dplyr::mutate(
      Virus = dplyr::recode(Virus,
                            EV71 = "EV-A71",
                            CVA6 = "CVA6",
                            EV68 = "EV-D68"
      )
    )
}

#### Misc. ####
make_phi_titer_summary <- function(raw_df, log_phi_draws, virus) {
  
  # use R-M method to calculate endpoint titer (TCID50), with boundary fixes
  assign_reed_muench_titer <- function(raw_df) {
    raw_df %>%
      dplyr::mutate(
        serumID = as.character(serumID),
        n_replicates = 2L
      ) %>%
      dplyr::group_by(serumID) %>%
      dplyr::group_modify(~{
        df <- dplyr::arrange(.x, dilutions)
        
        # lowest dilution boundary cases (8)
        if (df$outcome[1] == 0) {
          return(tibble::tibble(rm_titer = df$dilutions[1] / 2))
        }
        if (df$outcome[1] == 1) {
          return(tibble::tibble(rm_titer = df$dilutions[1]))
        }
        
        # highest dilution boundary cases (1024)
        if (df$outcome[nrow(df)] == 2) {
          return(tibble::tibble(rm_titer = df$dilutions[nrow(df)] * 2))
        }
        if (df$outcome[nrow(df)] == 1) {
          return(tibble::tibble(rm_titer = df$dilutions[nrow(df)]))
        }
        
        # standard Reed-Muench endpoint
        tibble::tibble(rm_titer = reed_muench_endpoint(df))
      }) %>%
      dplyr::ungroup()
  }
  
  rm_titers <- assign_reed_muench_titer(raw_df)
  
  # one row per serumID with its computed endpoint titre
  raw_one <- raw_df %>%
    dplyr::mutate(serumID = as.character(serumID)) %>%
    dplyr::distinct(serumID) %>%
    dplyr::arrange(serumID) %>%
    dplyr::left_join(rm_titers, by = "serumID")
  
  # posterior phi draws
  phi_mat <- exp(log_phi_draws)  # [draw x individual]
  
  stopifnot(nrow(raw_one) == ncol(phi_mat))
  
  # group by computed RM titre and pool ALL posterior draws across matching serumIDs
  out <- lapply(sort(unique(raw_one$rm_titer)), function(tt) {
    
    idx <- which(raw_one$rm_titer == tt)
    vals <- as.vector(phi_mat[, idx, drop = FALSE])
    
    dplyr::tibble(
      virus = virus,
      titer = tt,
      phi_med = stats::median(vals),
      phi_lo  = stats::quantile(vals, 0.025),
      phi_hi  = stats::quantile(vals, 0.975)
    )
  })
  
  dplyr::bind_rows(out) %>%
    dplyr::filter(is.finite(titer), titer >= 0)
}

make_phi_age_group_data <- function(serum_meta, log_phi_draws, virus) {
  stopifnot(nrow(serum_meta) == ncol(log_phi_draws))
  
  age_breaks <- c(0, 5, 10, 20, 35, 55, 80, 95)
  age_labels <- c("[0,5)", "[5,10)", "[10,20)", "[20,35)", "[35,55)", "[55,80)", "[80,95]")
  
  dplyr::tibble(
    age = as.numeric(serum_meta$age),
    phi = exp(matrixStats::colMedians(log_phi_draws)),
    virus = virus
  ) %>%
    dplyr::mutate(
      age_group = cut(
        age,
        breaks = age_breaks,
        labels = age_labels,
        right = FALSE,
        include.lowest = TRUE
      ),
      age_group = as.character(age_group),
      age_group = dplyr::if_else(age >= 80 & age <= 95, "[80,95]", age_group),
      age_group = factor(age_group, levels = age_labels)
    ) %>%
    dplyr::filter(!is.na(age_group), is.finite(phi), phi > 0)
}

make_conc_dilution_titer_data <- function(raw_df, log_phi_draws, virus) {
  raw_df <- raw_df %>%
    dplyr::mutate(serumID = as.character(serumID))
  
  serum_phi <- raw_df %>%
    dplyr::distinct(serumID) %>%
    dplyr::arrange(serumID)
  
  stopifnot(nrow(serum_phi) == ncol(log_phi_draws))
  
  phi_med <- exp(matrixStats::colMedians(log_phi_draws))
  
  serum_phi <- serum_phi %>%
    dplyr::mutate(phi = phi_med)
  
  raw_df %>%
    dplyr::left_join(serum_phi, by = "serumID") %>%
    dplyr::mutate(
      virus = virus,
      titer = suppressWarnings(as.numeric(titer)),
      dilutions = suppressWarnings(as.numeric(dilutions)),
      phi_over_d = phi / dilutions
    ) %>%
    dplyr::filter(
      is.finite(phi),
      is.finite(phi_over_d),
      is.finite(dilutions),
      is.finite(titer),
      phi > 0,
      phi_over_d > 0,
      dilutions > 0,
      titer > 0
    )
}