# R/01_antibody_ppc.R

library(dplyr)
library(tidyr)
library(stringr)
library(posterior)
library(ggplot2)
library(matrixStats)


# build survival archetypes from observed data
build_archetypes_from_observed <- function(raw_df, serum_id_col = "serumID",
                                           dilution_col = "dilutions",
                                           outcome_col = "outcome",
                                           n_replicates = NULL) {
  df <- raw_df %>%
    rename(
      serumID = all_of(serum_id_col),
      dilutions = all_of(dilution_col),
      outcome = all_of(outcome_col)
    ) %>%
    mutate(serumID = as.character(serumID))
  
  # infer n_replicates if not provided
  if (is.null(n_replicates)) {
    n_rep_cols <- length(grep("^rep", names(df)))
    n_replicates <- if (n_rep_cols > 0) n_rep_cols else max(df$outcome, na.rm = TRUE)
  }
  
  dat <- df %>%
    arrange(serumID, dilutions) %>%
    mutate(
      individual = as.integer(factor(serumID)),
      survival   = outcome / n_replicates
    )
  
  stopifnot(all(is.finite(dat$survival)))
  stopifnot(all(dat$survival >= 0 & dat$survival <= 1))
  
  d_seq <- sort(unique(dat$dilutions))
  n_individuals <- length(unique(dat$individual))
  
  patterns <- dat %>%
    arrange(individual, dilutions) %>%
    group_by(individual) %>%
    summarise(archetype = paste(survival, collapse = "-"), .groups = "drop")
  
  archetype_map <- patterns %>%
    distinct(archetype) %>%
    arrange(archetype) %>%
    mutate(archetype_id = row_number())
  
  patterns <- patterns %>% left_join(archetype_map, by = "archetype")
  
  archetype_counts <- patterns %>%
    dplyr::count(archetype_id, name = "n_individuals")
  
  ind_to_arch <- patterns %>% select(individual, archetype_id)
  
  obs_archetype <- dat %>%
    left_join(ind_to_arch, by = "individual") %>%
    group_by(archetype_id, dilutions) %>%
    summarise(mean_surv = mean(survival), .groups = "drop")
  
  facet_labels <- archetype_counts %>%
    arrange(archetype_id) %>%
    mutate(facet_label = paste0("Archetype ", archetype_id, "\n(n = ", n_individuals, ")"))
  
  list(
    dat = dat,
    d_seq = d_seq,
    n_individuals = n_individuals,
    n_replicates = n_replicates,
    patterns = patterns,
    ind_to_arch = ind_to_arch,
    obs_archetype = obs_archetype,
    facet_labels = facet_labels
  )
}

# extract posterior draws from antibody fits
extract_antibody_draws <- function(draws_obj, n_individuals, include_k0 = FALSE, thin_draws = NULL, seed = 1) {
  draws <- posterior::as_draws_df(draws_obj)
  
  if (!("k1" %in% names(draws))) stop("Expected parameter 'k1' in antibody draws.")
  
  if (include_k0) {
    if (!("k0" %in% names(draws))) stop("include_k0=TRUE but no 'k0' in draws.")
    k0 <- draws$k0
  } else {
    k0 <- rep(0, nrow(draws))
  }
  
  # try log_phi first
  logphi_cols <- paste0("log_phi[", 1:n_individuals, "]")
  phi_cols    <- paste0("phi[", 1:n_individuals, "]")
  
  has_logphi <- all(logphi_cols %in% names(draws))
  has_phi    <- all(phi_cols %in% names(draws))
  
  if (!has_logphi && !has_phi) {
    stop("Could not find log_phi[1:n] or phi[1:n] in antibody draws.")
  }
  
  if (has_logphi) {
    log_phi <- as.matrix(draws[, logphi_cols, drop = FALSE])
  } else {
    phi <- as.matrix(draws[, phi_cols, drop = FALSE])
    log_phi <- log(phi)
  }
  
  k1 <- draws$k1
  
  if (!is.null(thin_draws) && thin_draws < nrow(draws)) {
    set.seed(seed)
    idx <- sample.int(nrow(draws), thin_draws)
    k0 <- k0[idx]
    k1 <- k1[idx]
    log_phi <- log_phi[idx, , drop = FALSE]
  }
  
  list(k0 = k0, k1 = k1, log_phi = log_phi)
}

# prior draws
r_half_cauchy <- function(n, scale) abs(rcauchy(n, location = 0, scale = scale))

sample_prior_draws <- function(n_draws, n_individuals,
                               prior_k1  = function(n) r_half_cauchy(n, scale = 1),
                               prior_phi = function(n) r_half_cauchy(n, scale = 500),
                               include_k0 = FALSE,
                               prior_k0 = function(n) rnorm(n, mean = 0, sd = 1)) {
  k1 <- prior_k1(n_draws)
  k0 <- if (include_k0) prior_k0(n_draws) else rep(0, n_draws)
  
  phi_mat <- matrix(prior_phi(n_draws * n_individuals), nrow = n_draws, ncol = n_individuals)
  log_phi <- log(phi_mat)
  
  list(k0 = k0, k1 = k1, log_phi = log_phi)
}


# predictive simulation by archetype
# - for each archetype and dilution:
#   simulate z_rep and compute mean survival across individuals per draw

ppc_by_archetype_from_draws <- function(d_seq, n_replicates, ind_to_arch, draws_obj) {
  k0 <- draws_obj$k0
  k1 <- draws_obj$k1
  log_phi <- draws_obj$log_phi
  n_draws <- length(k1)
  
  out <- lapply(sort(unique(ind_to_arch$archetype_id)), function(aid) {
    inds <- ind_to_arch$individual[ind_to_arch$archetype_id == aid]
    
    surv_draws_by_dil <- sapply(d_seq, function(dval) {
      # p_mat is [draw x individual]
      p_mat <- pnorm(k0 + k1 * (log_phi[, inds, drop = FALSE] - log(dval)))
      
      z_rep <- matrix(
        rbinom(n = length(p_mat), size = n_replicates, prob = as.vector(p_mat)),
        nrow = n_draws
      )
      
      # mean survival across individuals per draw
      rowMeans(z_rep / n_replicates)
    })
    
    tibble(
      archetype_id = aid,
      dilution = d_seq,
      lo  = apply(surv_draws_by_dil, 2, quantile, 0.025),
      med = apply(surv_draws_by_dil, 2, quantile, 0.5),
      hi  = apply(surv_draws_by_dil, 2, quantile, 0.975)
    )
  })
  
  bind_rows(out)
}


# prior + posterior PPC overlay plot
antibody_ppc_overlay_plot <- function(virus_prefix,
                                      raw_df,
                                      draws_obj,
                                      thin_draws_post = 1000,
                                      ncol = 4,
                                      include_k0 = FALSE,
                                      
                                      # prior-line controls
                                      add_prior_lines = TRUE,
                                      n_prior_lines = 300,
                                      prior_alpha = 0.12,
                                      prior_linewidth = 0.25,
                                      prior_seed = 1,
                                      
                                      # priors
                                      prior_k1  = function(n) r_half_cauchy(n, scale = 1),
                                      prior_phi = function(n) r_half_cauchy(n, scale = 500),
                                      prior_k0  = function(n) rnorm(n, 0, 1),
                                      
                                      seed = 1) {
  
  arch <- build_archetypes_from_observed(raw_df, n_replicates = NULL)
  d_seq <- arch$d_seq
  n_individuals <- arch$n_individuals
  n_replicates <- arch$n_replicates
  
  # observed points
  obs <- arch$obs_archetype %>%
    dplyr::left_join(arch$facet_labels, by = "archetype_id") %>%
    dplyr::mutate(facet_label = factor(facet_label, levels = arch$facet_labels$facet_label))
  
  # posterior predictive
  post_draws <- extract_antibody_draws(
    draws_obj = draws_obj,
    n_individuals = n_individuals,
    include_k0 = include_k0,
    thin_draws = thin_draws_post,
    seed = seed
  )
  
  pp_post <- ppc_by_archetype_from_draws(d_seq, n_replicates, arch$ind_to_arch, post_draws) %>%
    dplyr::left_join(arch$facet_labels, by = "archetype_id") %>%
    dplyr::mutate(facet_label = factor(facet_label, levels = arch$facet_labels$facet_label))
  
  d_grid <- exp(seq(log(min(d_seq)), log(max(d_seq)), length.out = 120))
  
  prior_df <- simulate_prior_survival_lines(
    d_seq = d_grid,
    n_prior_samples = n_prior_lines,   
    k1_scale = 1,
    phi_scale = 500,
    seed = prior_seed
  )
  
  # replicate the same prior lines into every facet
  prior_facet_df <- tidyr::crossing(
    prior_df,
    facet_label = levels(pp_post$facet_label)
  ) %>%
    dplyr::mutate(
      facet_label = factor(facet_label, levels = levels(pp_post$facet_label))
    )
  
  ggplot() +
    # faint prior lines in the background
    { if (add_prior_lines)
      geom_line(
        data = prior_facet_df,
        aes(x = dilution, y = survival, group = draw),
        color = "grey60",
        alpha = prior_alpha,
        linewidth = prior_linewidth
      )
      else NULL
    } +
    # posterior predictive ribbon + median
    geom_ribbon(
      data = pp_post,
      aes(x = dilution, ymin = lo, ymax = hi),
      fill = "steelblue", alpha = 0.30
    ) +
    geom_line(
      data = pp_post,
      aes(x = dilution, y = med),
      color = "steelblue4", linewidth = 0.9
    ) +
    # observed points
    geom_point(
      data = obs,
      aes(x = dilutions, y = mean_surv),
      size = 2, color = "black"
    ) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~ facet_label, ncol = ncol) +
    labs(
      x = "Dilution",
      y = "Mean survival proportion",
      title = paste0(virus_prefix, ": Prior trajectories + Posterior predictive PPCs by archetype")
    ) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))
}

simulate_prior_survival_lines <- function(
    d_seq,
    n_prior_samples = 500,
    k1_scale = 1,
    phi_scale = 500,
    seed = 1
) {
  set.seed(seed)
  
  k1_prior  <- abs(rcauchy(n_prior_samples, location = 0, scale = k1_scale))
  phi_prior <- abs(rcauchy(n_prior_samples, location = 0, scale = phi_scale))
  log_phi   <- log(pmax(phi_prior, 1e-12))
  
  prior_df <- tidyr::crossing(
    draw = seq_len(n_prior_samples),
    dilution = d_seq
  ) %>%
    dplyr::mutate(
      survival = pnorm(k1_prior[draw] * (log_phi[draw] - log(dilution)))
    ) %>%
    dplyr::select(draw, dilution, survival)
  
  prior_df
}