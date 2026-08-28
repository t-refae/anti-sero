#### adapt existing sero_stan_data list for revised SIS/SIIS Stan files ####

add_sero_switches <- function(stan_data,
                              estimate_omega = 1L,
                              omega_fixed    = 0,
                              omega_prior_sd = 0.5,
                              sigma_by_state = 1L,
                              foi_piecewise  = 0L,
                              bracket_cuts   = c(0, 5, 10, 20, 40)) {
  
  need <- c("n", "age_max", "ages", "n_phi_draws", "log_phi_draws")
  stopifnot(all(need %in% names(stan_data)))
  
  if (missing(sigma_by_state) && !is.null(stan_data$sigma_by_state)) {
    sigma_by_state <- stan_data$sigma_by_state
  }
  
  n  <- stan_data$n
  M  <- stan_data$n_phi_draws
  lp <- stan_data$log_phi_draws
  
  if (n == M) {
    warning("n equals n_phi_draws; cannot infer orientation of log_phi_draws. ",
            "Check that individuals are in rows.")
  } else if (nrow(lp) == M && ncol(lp) == n) {
    lp <- t(lp)
  }
  
  stopifnot(nrow(lp) == n, ncol(lp) == M, !anyNA(lp))
  stopifnot(min(stan_data$ages) >= 0, max(stan_data$ages) <= stan_data$age_max)
  
  bracket_of_age <- if (foi_piecewise) {
    as.integer(cut(seq_len(stan_data$age_max),
                   breaks = c(bracket_cuts, Inf),
                   right = FALSE, labels = FALSE))
  } else {
    rep(0L, stan_data$age_max)
  }
  
  list(
    n              = as.integer(n),
    age_max        = as.integer(stan_data$age_max),
    ages           = as.array(as.integer(stan_data$ages)),
    n_phi_draws    = as.integer(M),
    log_phi_draws  = lp,
    estimate_omega = as.integer(estimate_omega),
    omega_fixed    = as.numeric(omega_fixed),
    omega_prior_sd = as.numeric(omega_prior_sd),
    sigma_by_state = as.integer(sigma_by_state),
    foi_piecewise  = as.integer(foi_piecewise),
    n_brackets     = if (foi_piecewise) length(bracket_cuts) else 0L,
    bracket_of_age = as.array(bracket_of_age),
    grainsize      = 1L
  )
}

#### model comparison: SI/SII/SIS/SIIS ####

find_log_lik_var <- function(draws, candidates = c("log_lik", "log_likelihood")) {
  vars <- posterior::variables(draws)
  base <- sub("\\[.*$", "", vars)
  present <- candidates[candidates %in% base]
  
  if (length(present) == 0L) {
    stop("None of the candidate log-likelihood variables (",
         paste(candidates, collapse = ", "), ") were found in this fit.",
         call. = FALSE)
  }
  if (length(present) > 1L) {
    stop("Fit contains more than one candidate log-likelihood variable (",
         paste(present, collapse = ", "),
         "). Name the intended one explicitly via log_lik_var.",
         call. = FALSE)
  }
  present
}

extract_log_lik_array <- function(draws, var = c("log_lik", "log_likelihood")) {
  d <- posterior::as_draws_array(draws)
  v <- find_log_lik_var(d, candidates = var)
  
  vars <- posterior::variables(d)
  sel  <- vars[sub("\\[.*$", "", vars) == v]   # original storage order preserved
  d    <- posterior::subset_draws(d, variable = sel)
  
  list(
    array = array(as.numeric(d), dim = dim(d)),
    var   = v
  )
}

# PSIS-LOO for a single fit
loo_from_draws <- function(draws, var = c("log_lik", "log_likelihood")) {
  ll    <- extract_log_lik_array(draws, var = var)
  r_eff <- loo::relative_eff(exp(ll$array))
  out   <- loo::loo(ll$array, r_eff = r_eff)
  attr(out, "log_lik_var") <- ll$var
  out
}

compare_sero_models_loo <- function(draws_list,
                                    virus,
                                    log_lik_var = c("log_lik", "log_likelihood"),
                                    k_threshold = 0.7) {
  
  if (is.null(names(draws_list)) || any(!nzchar(names(draws_list)))) {
    stop("draws_list must be a fully named list.", call. = FALSE)
  }
  if (length(draws_list) < 2L) {
    stop("Need at least two models to compare.", call. = FALSE)
  }
  
  pinned <- !is.null(names(log_lik_var)) && all(nzchar(names(log_lik_var)))
  if (pinned) {
    unpinned <- setdiff(names(draws_list), names(log_lik_var))
    if (length(unpinned)) {
      stop("log_lik_var is named but gives no entry for: ",
           paste(unpinned, collapse = ", "), call. = FALSE)
    }
  }
  
  loos <- lapply(names(draws_list), function(m) {
    v <- if (pinned) log_lik_var[[m]] else log_lik_var
    loo_from_draws(draws_list[[m]], var = v)
  })
  names(loos) <- names(draws_list)
  
  ## ensure consistent data for valid loo comparison
  n_obs <- vapply(loos, function(x) nrow(x$pointwise), integer(1))
  if (length(unique(n_obs)) != 1L) {
    stop(
      "Models have different numbers of pointwise observations (",
      paste(sprintf("%s=%d", names(n_obs), n_obs), collapse = ", "),
      "). LOO comparison is not valid.",
      call. = FALSE
    )
  }
  
  cmp <- loo::loo_compare(loos)
  cmp_df <- data.frame(
    model     = rownames(cmp),
    elpd_diff = as.numeric(cmp[, "elpd_diff"]),
    se_diff   = as.numeric(cmp[, "se_diff"]),
    stringsAsFactors = FALSE
  )
  
  est_df <- do.call(rbind, lapply(names(loos), function(m) {
    l <- loos[[m]]
    k <- l$diagnostics$pareto_k
    data.frame(
      model         = m,
      log_lik_var   = attr(l, "log_lik_var"),
      elpd_loo      = l$estimates["elpd_loo", "Estimate"],
      se_elpd_loo   = l$estimates["elpd_loo", "SE"],
      p_loo         = l$estimates["p_loo",    "Estimate"],
      looic         = l$estimates["looic",    "Estimate"],
      n_obs         = nrow(l$pointwise),
      max_pareto_k  = max(k),
      n_k_above_thr = sum(k > k_threshold),
      stringsAsFactors = FALSE
    )
  }))
  
  out <- merge(cmp_df, est_df, by = "model", sort = FALSE)
  out <- out[order(-out$elpd_loo), ]
  
  out$z_diff <- ifelse(out$se_diff > 0, out$elpd_diff / out$se_diff, NA_real_)
  
  out$virus <- virus
  out$best  <- out$model == out$model[1]
  
  rownames(out) <- NULL
  out[, c("virus", "model", "elpd_loo", "se_elpd_loo", "p_loo", "looic",
          "elpd_diff", "se_diff", "z_diff",
          "max_pareto_k", "n_k_above_thr", "n_obs", "log_lik_var", "best")]
}

combine_model_loo_tables <- function(...) {
  df <- dplyr::bind_rows(...)
  df$model <- factor(df$model, levels = c("SI", "SII", "SIS", "SIIS"))
  dplyr::arrange(df, virus, dplyr::desc(elpd_loo))
}



#### k0 sensitivity analysis ####
##
## refits antibody_mech_k0.stan with k0 held at -2, -1, 0, 1, 2
##
## panel  image
##
##   A    fitted dose-response curves
##   B    log(phi)
##   C    FOI posterior
##


## eta = k0 + k1 * (log_phi - log(d)).
k0_linpred <- function(k0, k1, log_phi, log_dil) {
  k0 + k1 * (log_phi - log_dil)
}

k0_label <- function(spec){
  if (identical(spec, "free")){
    "free" 
  } else {
    paste("k0 =", spec)
  }
}

K0_FOI_REGEX    <- "^lambda"


k0_compile <- function(stan_file = "Stan/antibody_mech_k0.stan",
                       dir = "Stan/bin") {
  ensure_dir(dir)
  cmdstanr::cmdstan_model(stan_file, dir = dir, cpp_options = list(stan_threads = TRUE))
}

fit_k0_stage1 <- function(model, stan_data, spec,
                          keep_log_phi_draws = FALSE, thin_to = NULL,
                          seed = 1, iter_warmup = 500, iter_sampling = 500,
                          chains = 4, parallel_chains = 4, refresh = 200,
                          cpp_options = stan_cpp, threads_per_chain = ab_threads) {
  
  free <- identical(spec, "free")
  zeta <- if (free) 0 else as.numeric(spec)
  
  stan_data$estimate_k0     <- as.integer(free)
  stan_data$transform_prior <- as.integer(!free)
  stan_data$k0_fixed        <- zeta
  stan_data$k0_prior_sd     <- K0_PRIOR_SD
  stan_data$phi_prior_scale <- K0_PHI_PRIOR_SCALE
  stan_data$k1_prior_scale  <- K1_PRIOR_SCALE
  
  fit <- model$sample(
    data = stan_data, seed = seed,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    chains = chains, parallel_chains = parallel_chains, refresh = refresh,
    # cpp_options = stan_cpp,
    threads_per_chain = ab_threads
  )
  
  dd <- posterior::as_draws_df(fit$draws())
  lp <- as.matrix(dplyr::select(dd, dplyr::starts_with("log_phi[")))
  ld <- as.matrix(dplyr::select(dd, dplyr::starts_with("log_LD50[")))
  
  out <- list(
    spec  = spec,
    zeta  = zeta,
    free  = free,
    label = k0_label(spec),
    d     = stan_data$d,
    scalars = tibble::tibble(
      k0           = dd$k0,
      k1           = dd$k1,
      mean_log_phi = matrixStats::rowMeans2(lp)
    ),
    log_phi_med = matrixStats::colMedians(lp),
    ld50_med    = matrixStats::colMedians(ld),
    diag = fit$summary(c("k0", "k1"), "rhat", "ess_bulk"),
    elpd = tryCatch(fit$loo("log_lik"), error = function(e) NULL)
  )
  
  if (keep_log_phi_draws) {
    out$log_phi_draws <- extract_log_phi_draws(
      draws_obj = dd, n_individuals = stan_data$n_individuals,
      thin_to = thin_to, seed = seed
    )
  }
  out
}

# returns FOI draws only for now
fit_k0_stage2 <- function(sero_stan_file, sero_stan_data, spec,
                          dir = "Stan/bin", seed = sero_seed,
                          iter_warmup = K0_S2_WARMUP,
                          iter_sampling = K0_S2_SAMPLING,
                          chains = K0_S2_CHAINS, parallel_chains = K0_S2_CHAINS,
                          refresh = K0_S2_REFRESH, 
                          # stan_cpp=stan_cpp,
                          threads_per_chain=sero_threads) {
  
  ensure_dir(dir)
  model <- cmdstanr::cmdstan_model(sero_stan_file, dir = dir, cpp_options = list(stan_threads = TRUE))
  
  fit <- model$sample(
    data = sero_stan_data, seed = seed,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    chains = chains, parallel_chains = parallel_chains, refresh = refresh,
    # stan_cpp=stan_cpp,
    threads_per_chain=sero_threads
  )
  
  vars <- grep(K0_FOI_REGEX, fit$metadata()$stan_variables, value = TRUE)
  stopifnot(length(vars) > 0)
  
  posterior::as_draws_df(fit$draws(variables = vars)) |>
    dplyr::select(-dplyr::any_of(c(".chain", ".iteration", ".draw"))) |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "param", values_to = "foi") |>
    dplyr::mutate(spec = spec, label = k0_label(spec))
}


k0_scale <- function() {
  ggplot2::scale_colour_viridis_d(end = 0.9, option = "C")
}

# fitted dose-response curves
plot_k0_curves_simple <- function(fits, n_serum = 3, n_grid = 60) {
  base <- fits[["0"]]
  ord  <- order(base$log_phi_med)
  sel  <- ord[round(seq(1, length(ord), length.out = n_serum))]
  grid <- seq(log(min(base$d)), log(max(base$d)), length.out = n_grid)
  
  d <- purrr::map_dfr(fits, function(f) {
    tidyr::expand_grid(serum = sel, log_dil = grid) |>
      dplyr::mutate(
        label   = f$label,
        log_phi = f$log_phi_med[serum],
        pi = stats::pnorm(k0_linpred(stats::median(f$scalars$k0),
                                     stats::median(f$scalars$k1),
                                     log_phi, log_dil))
      )
  })
  
  ggplot2::ggplot(d, ggplot2::aes(exp(log_dil), pi, colour = label)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(~ paste("serum", serum)) +
    ggplot2::scale_x_log10() +
    k0_scale() +
    ggplot2::labs(x = "dilution factor d",
                  y = expression(Phi(eta[i](d))), colour = NULL) +
    ggplot2::theme_bw(base_size = 9)
}

k0_max_curve_diff <- function(fits, n_grid = 200) {
  base <- fits[["0"]]
  grid <- seq(log(min(base$d)), log(max(base$d)), length.out = n_grid)
  
  purrr::map_dfr(fits, function(f) {
    tidyr::expand_grid(serum = seq_along(f$log_phi_med), log_dil = grid) |>
      dplyr::mutate(
        label = f$label,
        pi = stats::pnorm(k0_linpred(stats::median(f$scalars$k0),
                                     stats::median(f$scalars$k1),
                                     f$log_phi_med[serum], log_dil))
      )
  }) |>
    dplyr::group_by(serum, log_dil) |>
    dplyr::summarise(r = diff(range(pi)), .groups = "drop") |>
    dplyr::pull(r) |> max()
}

# log(phi) against k0
plot_k0_logphi <- function(fits) {
  d <- purrr::map_dfr(fits, function(f) {
    tibble::tibble(
      label   = f$label,
      zeta    = f$zeta,
      log_phi = f$log_phi_med
    )
  }) |>
    dplyr::mutate(
      label = factor(
        label,
        levels = c("k0 = -2", "k0 = -1", "k0 = 0", "k0 = 1", "k0 = 2")
      )
    )
  
  k1   <- stats::median(fits[["0"]]$scalars$k1)
  m0   <- stats::median(fits[["0"]]$log_phi_med)
  pred <- tibble::tibble(zeta = sort(unique(d$zeta))) |>
    dplyr::mutate(med = m0 - zeta / k1)
  
  ggplot2::ggplot(d, ggplot2::aes(factor(zeta), log_phi, colour = label)) +
    ggplot2::geom_boxplot(
      linewidth = 0.4,
      outlier.size = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = pred,
      ggplot2::aes(factor(zeta), med, group = 1),
      inherit.aes = FALSE,
      linetype = "dashed"
    ) +
    k0_scale() +
    ggplot2::labs(
      x = expression(k[0]),
      y = expression(median ~ log ~ phi[i])
    ) +
    ggplot2::theme_bw(base_size = 9)
}

# FOI posteriors
plot_k0_foi_simple <- function(foi_long) {
  
  foi_summary <- foi_long |>
    dplyr::filter(
      stringr::str_detect(param, "^lambda_long\\[\\d+\\]$")
    ) |>
    dplyr::mutate(
      param_index = as.integer(
        stringr::str_extract(param, "(?<=\\[)\\d+(?=\\])")
      )
    ) |>
    dplyr::group_by(spec, label, param_index) |>
    dplyr::summarise(
      median = stats::median(foi, na.rm = TRUE),
      lower  = unname(stats::quantile(foi, 0.025, na.rm = TRUE)),
      upper  = unname(stats::quantile(foi, 0.975, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(spec, param_index)
  
  foi_summary <- foi_summary |>
    dplyr::mutate(
      label = factor(
        label,
        levels = c("k0 = -2", "k0 = -1", "k0 = 0", "k0 = 1", "k0 = 2")
      )
    )
  
  ggplot2::ggplot(
    foi_summary,
    ggplot2::aes(
      x = param_index,
      group = label,
      colour = label
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = lower,
        ymax = upper,
        fill = ggplot2::after_scale(
          scales::alpha(colour, 0.3)
        )
      )#,
      # colour = NA,
      # show.legend = FALSE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = median,
        linetype = label
      ),
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      data = foi_summary |>
        dplyr::filter(param_index %% 5 == 1),
      ggplot2::aes(
        y = median,
        shape = label
      ),
      size = 1.4,
      stroke = 0.4
    ) +
    k0_scale() +
    ggplot2::scale_linetype_manual(
      values = c("solid", "longdash", "dashed", "dotdash", "dotted")
    ) +
    ggplot2::scale_shape_manual(
      values = c(16, 17, 15, 18, 3)
    ) +
    ggplot2::labs(
      x = "Age (Years)",
      y = "Force of Infection",
      colour = NULL,
      linetype = NULL,
      shape = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      legend.position = "bottom"
    ) +
    ggplot2::scale_x_log10()
}

k0_collect <- function(...) {
  fits <- list(...)
  stats::setNames(fits, vapply(fits, function(f) f$spec, character(1)))[K0_SPECS]
}

make_k0_figure_cva6 <- function(fits, foi_long, specs = K0_FIXED) {
  f <- fits[specs]
  
  p_a <- plot_k0_curves_simple(f) +
    ggplot2::guides(colour = "none")
  
  p_b <- plot_k0_logphi(f) +
    ggplot2::guides(colour = "none")
  
  p_c <- plot_k0_foi_simple(foi_long)
  
  p <- p_a / (p_b | p_c)
  
  p +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(legend.position = "bottom")
}

make_k0_table <- function(fits, virus) {
  purrr::map_dfr(fits, function(f) {
    tibble::tibble(
      virus         = virus,
      k0_spec       = f$spec,
      k0            = sprintf("%.3f", stats::median(f$scalars$k0)),
      k1            = sprintf("%.3f", stats::median(f$scalars$k1)),
      median_log_phi  = sprintf("%.3f", stats::median(f$log_phi_med)),
      median_log_LD50 = sprintf("%.3f", stats::median(f$ld50_med)),
      elpd_loo = if (is.null(f$elpd)) NA_character_ else
        sprintf("%.2f (%.2f)",
                f$elpd$estimates["elpd_loo", "Estimate"],
                f$elpd$estimates["elpd_loo", "SE"]),
      max_rhat = sprintf("%.3f", max(f$diag$rhat)),
      min_ess  = sprintf("%.0f", min(f$diag$ess_bulk))
    )
  })
}

# for free estimation of k0, include later if nec.
plot_k0_ridge <- function(fit_free, virus = "CVA6") {
  d   <- fit_free$scalars
  k1  <- stats::median(d$k1)
  a   <- stats::median(d$mean_log_phi) + stats::median(d$k0) / k1
  
  ggplot2::ggplot(d, ggplot2::aes(k0, mean_log_phi)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.5) +
    ggplot2::geom_abline(intercept = a, slope = -1 / k1,
                         linewidth = 0.7, linetype = "dashed") +
    ggplot2::labs(
      x = expression(k[0]), y = expression(bar(log~phi)),
      subtitle = sprintf(
        "%s | correlation %.3f | posterior SD %.2f vs prior SD %.0f | min ESS %.0f",
        virus, stats::cor(d$k0, d$mean_log_phi), stats::sd(d$k0),
        K0_PRIOR_SD, min(fit_free$diag$ess_bulk))
    ) +
    ggplot2::theme_bw(base_size = 9)
}



#### mean vs median phi (clarifying plots) ####

mixture_phi_median <- function(w, mu_x, sigma) {
  k  <- length(mu_x)
  sg <- if (length(sigma) == 1L) rep(sigma, k) else sigma[seq_len(k)]
  w  <- w[seq_len(k)]
  if (sum(w) <= 0 || anyNA(w)) return(NA_real_)
  w <- w / sum(w)
  f <- function(z) sum(w * stats::pnorm((z - mu_x) / sg)) - 0.5
  exp(stats::uniroot(f, lower = min(mu_x) - 8 * max(sg),
                     upper = max(mu_x) + 8 * max(sg))$root)
}

mixture_phi_summaries <- function(probs_age, mu_x, sigma) {
  k  <- length(mu_x)
  sg <- if (length(sigma) == 1L) rep(sigma, k) else sigma[seq_len(k)]
  P  <- probs_age[, seq_len(k), drop = FALSE]
  P  <- P / rowSums(P)
  data.frame(
    Age    = seq_len(nrow(P)) - 1L,
    mean   = as.numeric(P %*% exp(mu_x + 0.5 * sg^2)),
    median = apply(P, 1, mixture_phi_median, mu_x = mu_x, sigma = sigma)
  )
}

#### seroreversion plots ####

# recreates figs 4/5/6 for SIS/SIIS fits

sr_model_label <- function(model) {
  switch(toupper(model),
         SIS  = "SIS",
         SIIS = "SIIS",
         toupper(model))
}

sr_title <- function(virus, model) {
  paste0(pretty_virus(virus), " (", sr_model_label(model), ")")
}

sr_get_draw <- function(draws_df, name) {
  cand <- c(name, paste0(name, "[1]"))
  hit  <- cand[cand %in% names(draws_df)]
  if (!length(hit)) return(NULL)
  as.numeric(draws_df[[hit[1]]])
}

sr_get_summary <- function(summary_df, name) {
  cand <- c(name, paste0(name, "[1]"))
  out  <- summary_df %>% filter(variable %in% cand) %>% pull(mean)
  if (!length(out)) NULL else out[1]
}

sr_n_states <- function(draws_df) {
  if ("mu[3]" %in% names(draws_df) || "mu_x[3]" %in% names(draws_df)) 3L else 2L
}

sr_state_prob_mats <- function(draws_df, age_max, n_states) {
  lapply(seq_len(n_states), function(j) {
    cols <- paste0("state_probs[", seq_len(age_max + 1L), ",", j, "]")
    miss <- setdiff(cols, names(draws_df))
    if (length(miss)) {
      stop("state_probs missing from draws (", length(miss), " columns, e.g. ",
           miss[1], "). Was the fit run with the revised SIS/SIIS Stan file?")
    }
    as.matrix(draws_df[, cols, drop = FALSE])
  }) %>% setNames(serostate_levels[seq_len(n_states)])
}

sr_mean_probs <- function(prob_mats) {
  P <- vapply(prob_mats, colMeans, numeric(ncol(prob_mats[[1]])))
  if (ncol(P) == 2L) P <- cbind(P, Ip = 0)
  colnames(P) <- serostate_levels
  P
}

sr_mixture_summaries <- function(probs_age, mu_x, sigma) {
  k  <- length(mu_x)
  sg <- if (length(sigma) == 1L) rep(sigma, k) else sigma[seq_len(k)]
  P  <- probs_age[, seq_len(k), drop = FALSE]
  P  <- P / rowSums(P)
  
  med <- apply(P, 1, function(w) {
    if (any(is.na(w))) return(NA_real_)
    f <- function(z) sum(w * stats::pnorm((z - mu_x) / sg)) - 0.5
    exp(stats::uniroot(f,
                       lower = min(mu_x) - 8 * max(sg),
                       upper = max(mu_x) + 8 * max(sg))$root)
  })
  
  data.frame(
    Age    = seq_len(nrow(P)) - 1L,
    mean   = as.numeric(P %*% exp(mu_x + 0.5 * sg^2)),
    median = med
  )
}


## fig 4

plot_sr_weighted_avg_phi <- function(virus, model, plot_data,
                                     draws_df, summary_df,
                                     seed = 1, eps = 0.95) {
  
  n_states <- sr_n_states(draws_df)
  age_max  <- plot_data$age_max
  n        <- length(plot_data$phi_median)
  
  prob_mat <- extract_prob_matrix(summary_df, n = n, k = n_states)
  
  set.seed(seed)
  state <- integer(n); below <- logical(n)
  for (i in seq_len(n)) {
    p <- prob_mat[i, ]; p[is.na(p)] <- 0; p <- p / sum(p)
    state[i] <- sample.int(n_states, 1, prob = p)
    below[i] <- (p[state[i]] < eps)
  }
  
  pts <- data.frame(
    age   = plot_data$ages,
    phi   = plot_data$phi_median,
    state = factor(serostate_levels[state], levels = serostate_levels),
    conf  = ifelse(below, "<95%", ">95%")
  )
  
  mu_x <- summary_df %>%
    filter(grepl("^mu_x\\[", variable)) %>%
    arrange(variable) %>% pull(mean)
  
  sigma <- summary_df %>%
    filter(grepl("^sigma_phi\\[", variable) | variable == "sigma") %>%
    arrange(variable) %>% pull(mean)
  
  probs_age <- sr_mean_probs(sr_state_prob_mats(draws_df, age_max, n_states))
  
  line_df <- sr_mixture_summaries(probs_age, mu_x, sigma) %>%
    tidyr::pivot_longer(c(mean, median),
                        names_to = "Summary", values_to = "phi") %>%
    mutate(Summary = factor(Summary, c("mean", "median"), c("Mean", "Median")))
  
  ggplot() +
    geom_point(data = pts, aes(x = age, y = phi, color = state, shape = conf),
               alpha = 0.85, show.legend = TRUE) +
    geom_line(data = line_df, aes(x = Age, y = phi, linetype = Summary),
              linewidth = 1.2) +
    scale_linetype_manual(values = c(Mean = "solid", Median = "dashed")) +
    scale_y_log10(breaks = c(10, 30, 100, 300, 1000, 3000)) +
    scale_serostate_color(state_names = serostate_levels, drop = FALSE) +
    scale_shape_manual(values = c(">95%" = 16, "<95%" = 0)) +
    scale_x_continuous(breaks = seq(0, 90, by = 10), limits = c(0, age_max)) +
    labs(title = sr_title(virus, model), x = "Age (Years)", y = expression(phi),
         color = "Serostate", linetype = NULL, shape = "Confidence") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 36, angle = 0, vjust = 0.5, face = "bold"),
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14)
    ) +
    guides(shape = "none",
           color = guide_legend(override.aes = list(size = 7)),
           linetype = guide_legend(order = 2))
}


## fig 5

plot_sr_serodynamics <- function(virus, model, plot_data, draws_df) {
  
  n_states  <- sr_n_states(draws_df)
  age_max   <- plot_data$age_max
  prob_mats <- sr_state_prob_mats(draws_df, age_max, n_states)
  state_names <- serostate_levels[seq_len(n_states)]
  
  df_long <- purrr::imap_dfr(prob_mats, function(M, nm) {
    data.frame(
      Age   = seq_len(ncol(M)) - 1L,
      State = nm,
      mean  = colMeans(M),
      lo    = apply(M, 2, stats::quantile, 0.025),
      hi    = apply(M, 2, stats::quantile, 0.975)
    )
  }) %>%
    mutate(State = factor(State, levels = state_names))
  
  legend_df <- data.frame(Age = NA_real_, mean = NA_real_,
                          lo = NA_real_, hi = NA_real_,
                          State = factor(state_names, levels = state_names))
  
  ggplot(df_long, aes(x = Age)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = State), alpha = 0.18) +
    geom_line(aes(y = mean, color = State), linewidth = 1.3) +
    geom_ribbon(data = legend_df, aes(ymin = lo, ymax = hi, fill = State),
                alpha = 0.18, show.legend = TRUE, na.rm = TRUE) +
    geom_line(data = legend_df, aes(y = mean, color = State),
              linewidth = 1.3, show.legend = TRUE, na.rm = TRUE) +
    scale_serostate_color(state_names = state_names,
                          limits = state_names, drop = FALSE) +
    scale_serostate_fill(state_names = state_names,
                         limits = state_names, drop = FALSE) +
    scale_x_continuous(breaks = seq(0, 90, by = 10), limits = c(0, age_max)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = sr_title(virus, model), x = "Age (Years)", y = "Proportion",
         color = "Serostate", fill = "Serostate") +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      legend.box = "vertical", legend.direction = "vertical",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 14),
      legend.key.height = grid::unit(1, "cm"),
      legend.key.width  = grid::unit(2, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))
}


## fig 6

plot_sr_foi_new_old <- function(virus, model, plot_data, draws_df,
                                old_age_offset = 1) {
  
  age_max <- plot_data$age_max
  lam_cols <- paste0("lambda_long[", seq_len(age_max), "]")
  stopifnot(all(lam_cols %in% names(draws_df)))
  FOI_mat <- as.matrix(draws_df[, lam_cols, drop = FALSE])   # draws x age_max
  
  new_df <- data.frame(
    Age  = seq_len(age_max),
    mean = colMeans(FOI_mat),
    q5   = apply(FOI_mat, 2, stats::quantile, 0.05),
    q95  = apply(FOI_mat, 2, stats::quantile, 0.95)
  )
  
  old_fit  <- readRDS(file.path("data", "processed",
                                paste0(virus, "_model5_fit.rds")))
  old_post <- rstan::extract(old_fit, pars = c("lambda", "beta"))
  age_seq  <- seq_len(age_max)
  
  old_FOI_mat <- vapply(
    seq_along(old_post$beta),
    function(i) old_post$lambda[i] *
      exp(-old_post$beta[i] * (age_seq - old_age_offset)),
    numeric(length(age_seq))
  )
  
  old_df <- data.frame(
    Age  = age_seq,
    mean = rowMeans(old_FOI_mat),
    q5   = matrixStats::rowQuantiles(old_FOI_mat, probs = 0.05),
    q95  = matrixStats::rowQuantiles(old_FOI_mat, probs = 0.95)
  )
  
  virus_col <- virus_cols[[virus]]
  
  ggplot() +
    geom_ribbon(data = new_df, aes(x = Age, ymin = q5, ymax = q95, fill = "New"),
                alpha = 0.3) +
    geom_line(data = new_df, aes(x = Age, y = mean, color = "New"),
              linewidth = 1.5) +
    geom_ribbon(data = old_df, aes(x = Age, ymin = q5, ymax = q95, fill = "Old"),
                alpha = 0.3) +
    geom_line(data = old_df, aes(x = Age, y = mean, color = "Old"),
              linewidth = 1.5) +
    labs(title = sr_title(virus, model), x = "Age (Years)",
         y = "Force of Infection",
         color = "Force of Infection", fill = "Force of Infection") +
    scale_x_log10() +
    scale_y_continuous(breaks = seq(0, 0.6, by = 0.2), limits = c(0, 0.66)) +
    scale_color_manual(values = c(New = virus_col, Old = "grey80")) +
    scale_fill_manual(values  = c(New = virus_col, Old = "grey80")) +
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title   = element_text(size = 16, face = "bold"),
      axis.text    = element_text(size = 12),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text  = element_text(size = 16),
      legend.position = "right",
      legend.key.width = grid::unit(3, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.5)))
}


plot_sr_omega <- function(virus, model, draws_df, omega_prior_sd = 0.1,
                          n_prior = 4000, seed = 1) {
  
  om <- sr_get_draw(draws_df, "omega")
  if (is.null(om) || stats::sd(om) < 1e-12) {
    return(ggplot() + theme_void() +
             labs(title = paste0(sr_title(virus, model), ": omega fixed")))
  }
  
  set.seed(seed)
  d <- bind_rows(
    data.frame(omega = om, dist = "Posterior"),
    data.frame(omega = abs(stats::rnorm(n_prior, 0, omega_prior_sd)),
               dist = "Prior")
  )
  
  hl <- stats::median(om)
  hl <- if (hl > 0) log(2) / hl else NA_real_
  
  ggplot(d, aes(x = omega, fill = dist)) +
    annotate("rect", xmin = 0, xmax = log(2) / 20,
             ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "grey40") +
    geom_density(alpha = 0.45, colour = NA) +
    scale_fill_manual(values = c(Prior = "grey60", Posterior = virus_cols[[virus]])) +
    labs(
      title = sr_title(virus, model),
      subtitle = paste0(
        "Shaded: seroreversion half-life > 20 years.  ",
        "Posterior median half-life = ", sprintf("%.1f", hl), " years"),
      x = expression(omega~"(per year)"), y = "Density", fill = NULL
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "bottom")
}


## table S2

make_sr_param_table <- function(virus, model, summary_df) {
  
  n_states <- if (any(grepl("^mu_x\\[3\\]", summary_df$variable))) 3L else 2L
  
  recode_map <- c(
    "mu_x[1]"      = "\\mu_{S}",
    "mu_x[2]"      = "\\mu_{I_{++}}",
    "mu_x[3]"      = "\\mu_{I_{+}}",
    "sigma_phi[1]" = "\\sigma_{\\phi,S}",
    "sigma_phi[2]" = "\\sigma_{\\phi,I_{++}}",
    "sigma_phi[3]" = "\\sigma_{\\phi,I_{+}}",
    "psi"          = "\\psi",
    "omega"        = "\\omega",
    "lambda_1[1]"  = "\\lambda_{1}",
    "kappa[1]"     = "\\kappa"
  )
  
  vars <- c("mu_x[1]", "mu_x[2]", if (n_states == 3) "mu_x[3]",
            "sigma_phi[1]", "sigma_phi[2]", if (n_states == 3) "sigma_phi[3]",
            if (n_states == 3) "psi",
            "omega", "lambda_1[1]", "kappa[1]")
  vars <- vars[vars %in% summary_df$variable]
  
  tab <- summary_df %>%
    filter(variable %in% vars) %>%
    select(variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail) %>%
    mutate(variable = recode(variable, !!!recode_map))
  
  tab$variable <- factor(tab$variable, levels = recode_map[vars])
  tab <- tab %>% arrange(variable)
  tab$variable <- paste0("$", tab$variable, "$")
  
  gt::gt(tab) %>%
    gt::tab_header(title = paste0(sr_title(virus, model),
                                  ": Summary Statistics of Fitted Parameters")) %>%
    gt::fmt_number(columns = c(mean, median, sd, q5, q95), decimals = 3) %>%
    gt::fmt_number(columns = rhat, decimals = 3) %>%
    gt::fmt_number(columns = c(ess_bulk, ess_tail), decimals = 0) %>%
    gt::fmt_markdown(columns = variable) %>%
    gt::cols_label(variable = "VARIABLE", mean = "MEAN", median = "MEDIAN",
                   sd = "SD", q5 = "Q5", q95 = "Q95", rhat = "RHAT",
                   ess_bulk = "ESS_BULK", ess_tail = "ESS_TAIL")
}


sr_panel <- function(what, virus, model, plot_data, draws_df, summary_df,
                     omega_prior_sd = 0.1) {
  switch(what,
         weighted_avg_phi = plot_sr_weighted_avg_phi(virus, model, plot_data,
                                                     draws_df, summary_df),
         serodynamics     = plot_sr_serodynamics(virus, model, plot_data, draws_df),
         FOI_new_old      = plot_sr_foi_new_old(virus, model, plot_data, draws_df),
         omega            = plot_sr_omega(virus, model, draws_df, omega_prior_sd),
         stop("unknown panel type: ", what)
  )
}

save_sr_combined <- function(what, model,
                             pd_cva6, dr_cva6, sm_cva6,
                             pd_ev71, dr_ev71, sm_ev71,
                             pd_ev68, dr_ev68, sm_ev68,
                             dir = "outputs/seroreversion",
                             omega_prior_sd = 0.5) {
  
  p_cva6 <- sr_panel(what, "CVA6", model, pd_cva6, dr_cva6, sm_cva6, omega_prior_sd)
  p_ev71 <- sr_panel(what, "EV71", model, pd_ev71, dr_ev71, sm_ev71, omega_prior_sd)
  p_ev68 <- sr_panel(what, "EV68", model, pd_ev68, dr_ev68, sm_ev68, omega_prior_sd)
  
  grid <- switch(what,
                 FOI_new_old      = plot_patchwork(p_cva6, p_ev71, p_ev68, FOI = TRUE),
                 weighted_avg_phi = plot_patchwork(p_cva6, p_ev71, p_ev68),
                 plot_patchwork(p_cva6, p_ev71, p_ev68, lines = FALSE)
  )
  
  dims <- if (identical(what, "omega")) c(12, 8) else c(12, 10)
  
  ensure_dir(dir)
  out <- file.path(dir, paste0("combined_", tolower(sr_model_label(model)),
                               "_", what, ".pdf"))
  save_plot_pdf(grid, out, dims[1], dims[2])
  out
}

save_sr_param_tables <- function(specs,
                                 dir  = "outputs/seroreversion",
                                 file = "combined_seroreversion_parameter_tables.pdf") {
  ensure_dir(dir)
  tmp <- file.path(tempdir(), "sr_tables")
  ensure_dir(tmp)
  
  paths <- vapply(specs, function(s) {
    save_gt_pdf(
      make_sr_param_table(s$virus, s$model, s$summary),
      file.path(tmp, paste0(s$virus, "_", tolower(s$model), "_param_table.pdf"))
    )
  }, character(1))
  
  out <- file.path(dir, file)
  combine_pdfs(paths, out)
  unlink(paths)
  out
}


#### vectorised sigma (serostate-specific) ####

serostate_reclassification <- function(summary_a, summary_b, n, k) {
  modal <- function(s) apply(extract_prob_matrix(s, n = n, k = k), 1, which.max)
  a <- modal(summary_a); b <- modal(summary_b)
  list(prop_changed = mean(a != b, na.rm = TRUE),
       table = table(shared = a, by_state = b))
}

#### convergence diagnostics ####

# one row per fit
sero_convergence_row <- function(summary_df, draws_df, virus, model,
                                 sigma_spec = "by state", diagnostics = NULL) {
  
  key <- summary_df |>
    dplyr::filter(grepl(
      "^(mu|mu_x|sigma_phi|sigma_x|sig_x|psi|omega|lambda_1|kappa)(\\[|$)",
      variable))
  
  n_states <- if (any(grepl("^mu_x\\[3\\]", summary_df$variable))) 3L else 2L
  sg <- grep("^sigma_phi\\[", names(draws_df), value = TRUE)
  
  # per-chain minima detect a single stuck chain that pooled minima would hide
  per_chain <- draws_df |>
    dplyr::group_by(.chain) |>
    dplyr::summarise(ms = min(dplyr::across(dplyr::all_of(sg))),
                     .groups = "drop")
  
  om_free <- "omega" %in% names(draws_df) &&
    stats::sd(draws_df$omega, na.rm = TRUE) > 1e-10
  
  tibble::tibble(
    virus = virus, model = model, n_states = n_states, sigma_spec = sigma_spec,
    omega = if (om_free) "free" else "fixed",
    n_chains = length(unique(draws_df$.chain)),
    n_draws  = nrow(draws_df),
    max_rhat     = max(key$rhat, na.rm = TRUE),
    n_rhat_gt101 = sum(key$rhat > 1.01, na.rm = TRUE),
    min_ess_bulk = min(key$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(key$ess_tail, na.rm = TRUE),
    min_sigma          = min(per_chain$ms),
    n_chains_collapsed = sum(per_chain$ms < 0.2),
    lp_range = if ("lp__" %in% names(draws_df)) {
      m <- tapply(draws_df$lp__, draws_df$.chain, mean); max(m) - min(m)
    } else NA_real_,
    n_divergent = if (!is.null(diagnostics)) sum(diagnostics$divergent__) else NA_integer_,
    converged = max(key$rhat, na.rm = TRUE) < 1.01 &
      min(key$ess_bulk, na.rm = TRUE) > 400
  )
}

label_convergence_failures <- function(df) {
  df |>
    dplyr::mutate(diagnosis = dplyr::case_when(
      converged                                    ~ "converged",
      n_chains_collapsed > 0 & omega == "free"     ~ "sigma collapse + omega ridge",
      n_chains_collapsed > 0                       ~ "sigma collapse",
      omega == "free"                              ~ "omega ridge",
      TRUE                                         ~ "other")) |>
    dplyr::arrange(converged, n_states, virus, model)
}

make_convergence_table <- function(df) {
  gt::gt(df |> dplyr::select(virus, model, n_states, sigma_spec, omega,
                             max_rhat, min_ess_bulk, min_sigma,
                             n_chains_collapsed, diagnosis)) |>
    gt::tab_header(title = "Stage 2 convergence diagnostics",
                   subtitle = "Convergence: R-hat < 1.01 and bulk ESS > 400 on all model parameters") |>
    gt::fmt_number(columns = c(max_rhat, min_sigma), decimals = 3) |>
    gt::fmt_number(columns = min_ess_bulk, decimals = 0) |>
    gt::data_color(columns = diagnosis,
                   fn = function(x) ifelse(x == "converged", "#d9ead3", "#f4cccc")) |>
    gt::cols_label(n_states = "STATES", sigma_spec = "SIGMA", omega = "OMEGA",
                   max_rhat = "MAX RHAT", min_ess_bulk = "MIN ESS",
                   min_sigma = "MIN SIGMA", n_chains_collapsed = "CHAINS COLLAPSED")
}

#### stratify results by calendar year ####

# age x year summary
stationarity_descriptive <- function(serum_meta, log_phi_draws, virus,
                                     cuts = c(0,1,2,seq(5,95, by=5))) {
  tibble::tibble(
    age    = as.integer(round(serum_meta$age)),
    survey = as.integer(serum_meta$year),
    phi    = exp(matrixStats::colMedians(log_phi_draws))
  ) |>
    dplyr::mutate(band = cut(age, cuts, right = FALSE)) |>
    dplyr::group_by(band, survey) |>
    dplyr::summarise(n = dplyr::n(),
                     median_phi = median(phi),
                     mean_log_phi = mean(log(phi)),
                     se = sd(log(phi)) / sqrt(dplyr::n()),
                     .groups = "drop") |>
    dplyr::mutate(virus = virus)
}

plot_stationarity_descriptive <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(band, mean_log_phi,
                                   colour = factor(survey), group = survey)) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = mean_log_phi - 1.96 * se,
                                          ymax = mean_log_phi + 1.96 * se),
                             position = ggplot2::position_dodge(0.3)) +
    ggplot2::geom_line(position = ggplot2::position_dodge(0.3)) +
    ggplot2::facet_wrap(~virus, scales = "free_y") +
    ggplot2::labs(x = "Age band", y = expression(mean~log~phi), colour = "Survey") +
    ggplot2::theme_bw()
}

subset_sero_data_by_year <- function(stan_data, serum_meta, years) {
  keep <- as.integer(serum_meta$year) %in% years
  stopifnot(sum(keep) > 50)
  lp <- stan_data$log_phi_draws
  stan_data$log_phi_draws <- if (nrow(lp) == stan_data$n_phi_draws) {
    lp[, keep, drop = FALSE]                      # draws x individuals
  } else {
    lp[keep, , drop = FALSE]                      # individuals x draws
  }
  stan_data$ages <- stan_data$ages[keep]
  stan_data$n    <- sum(keep)
  if (!is.null(stan_data$grainsize)) stan_data$grainsize <- 1L
  stan_data
}

strata_contrast <- function(draws_a, draws_b, label_a, label_b,
                            pars = c("lambda_1", "kappa", "psi"),
                            n_draw = 4000, seed = 1) {
  set.seed(seed)
  pars <- intersect(pars, intersect(names(draws_a), names(draws_b)))
  purrr::map_dfr(pars, function(p) {
    a <- as.numeric(draws_a[[p]]); b <- as.numeric(draws_b[[p]])
    d <- sample(b, n_draw, TRUE) - sample(a, n_draw, TRUE)
    tibble::tibble(
      parameter   = p,
      median_a    = stats::median(a),
      median_b    = stats::median(b),
      ratio       = stats::median(b) / stats::median(a),
      diff_median = stats::median(d),
      diff_q5     = stats::quantile(d, 0.05),
      diff_q95    = stats::quantile(d, 0.95),
      p_increase  = mean(d > 0),
      stratum_a   = label_a,
      stratum_b   = label_b
    )
  })
}

strata_contrast_all <- function(draws_list,
                                pars = c("lambda_1", "kappa", "psi"),
                                n_draw = 4000, seed = 1) {
  labs <- names(draws_list)
  pairs <- c(
    lapply(seq_len(length(labs) - 1), function(i) c(i, i + 1)),   # adjacent
    if (length(labs) > 2) list(c(1, length(labs)))                 # first vs last
  )
  purrr::map_dfr(pairs, function(ij) {
    strata_contrast(draws_list[[ij[1]]], draws_list[[ij[2]]],
                    labs[ij[1]], labs[ij[2]], pars, n_draw, seed)
  })
}


strata_trend <- function(draws_list, par = "lambda_1", n_draw = 4000, seed = 1) {
  set.seed(seed)
  X <- vapply(draws_list, function(d) sample(as.numeric(d[[par]]), n_draw, TRUE),
              numeric(n_draw))
  tibble::tibble(
    parameter        = par,
    p_monotone_up    = mean(X[, 1] < X[, 2] & X[, 2] < X[, 3]),
    p_monotone_down  = mean(X[, 1] > X[, 2] & X[, 2] > X[, 3]),
    p_last_gt_first  = mean(X[, 3] > X[, 1]),
    ratio_last_first = stats::median(X[, 3]) / stats::median(X[, 1])
  )
}

strata_seroprev <- function(draws, age_max, label, ages = c(1, 5, 10, 20, 40)) {
  lam <- as.matrix(draws[, paste0("lambda_long[", seq_len(age_max), "]")])
  # S(a) = exp(-cumsum(lambda)) for the SI model
  S <- t(apply(lam, 1, function(x) exp(-cumsum(x))))
  purrr::map_dfr(ages, function(a) {
    p <- 1 - S[, a]
    tibble::tibble(age = a, stratum = label,
                   median = stats::median(p),
                   q5 = stats::quantile(p, 0.05),
                   q95 = stats::quantile(p, 0.95))
  })
}

plot_strata_foi <- function(virus, draws_list, age_max) {
  d <- purrr::imap_dfr(draws_list, function(dr, lab) {
    M <- as.matrix(dr[, paste0("lambda_long[", seq_len(age_max), "]")])
    tibble::tibble(Age = seq_len(age_max),
                   median = matrixStats::colMedians(M),
                   lo = matrixStats::colQuantiles(M, probs = 0.05),
                   hi = matrixStats::colQuantiles(M, probs = 0.95),
                   stratum = lab)
  })
  ggplot2::ggplot(d, ggplot2::aes(Age, median, colour = stratum, fill = stratum)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Age (Years)", y = "Force of Infection", title = virus,
                  colour = NULL, fill = NULL) +
    ggplot2::theme_minimal()
}


#### k1 / phi prior sensitivity analyses ####

PS_LEVELS <- c("baseline", "replicate", "phi tight", "phi wide", "k1 tight", "k1 wide")

ps_label_factor <- function(x) factor(x, levels = PS_LEVELS)

# censoring
ps_censoring_summary <- function(prep, virus) {
  z <- prep$stan_data$z
  R <- prep$stan_data$n_replicates
  M <- prep$stan_data$n_dilutions
  d <- prep$stan_data$d
  
  stopifnot(!is.unsorted(d))          # z[, M] must be the highest dilution
  right <- z[, M] == R                # still fully neutralising at the highest
  left  <- z[, 1] == 0                # no neutralisation even at the lowest
  
  tibble::tibble(
    virus         = virus,
    n             = nrow(z),
    d_min         = min(d),
    d_max         = max(d),
    n_right_cens  = sum(right),
    n_left_cens   = sum(left),
    n_censored    = sum(right | left),
    prop_censored = mean(right | left)
  )
}

# stage 1

ps_sensitivity_df <- function(x) {
  if (is.data.frame(x)) tibble::as_tibble(x) else tibble::as_tibble(x$sensitivity)
}

fit_ps_stage1 <- function(model, stan_data, label, phi_scale, k1_scale,
                          thin_to = NULL,
                          sens_vars = c("k1", "mean_log_phi"),
                          seed = 1, iter_warmup = 500, iter_sampling = 500,
                          chains = 4, parallel_chains = 4, refresh = 200,
                          threads_per_chain = 1) {
  
  stan_data$estimate_k0     <- 0L
  stan_data$k0_fixed        <- 0
  stan_data$k0_prior_sd     <- K0_PRIOR_SD
  stan_data$phi_prior_scale <- phi_scale
  stan_data$k1_prior_scale  <- k1_scale
  
  fit <- model$sample(
    data = stan_data, seed = seed,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    chains = chains, parallel_chains = parallel_chains, refresh = refresh,
    threads_per_chain = threads_per_chain
  )
  
  # priorsense: sensitivity of each quantity to each prior component
  sens <- purrr::map_dfr(
    list(phi = "phi", k1 = "k1", joint = NULL),
    function(sel) {
      ps_sensitivity_df(
        priorsense::powerscale_sensitivity(
          fit, variable = sens_vars, prior_selection = sel
        )
      )
    },
    .id = "component"
  ) |>
    dplyr::mutate(label = label, .before = 1)
  
  dd <- posterior::as_draws_df(fit$draws())
  lp <- as.matrix(dplyr::select(dd, dplyr::starts_with("log_phi[")))
  
  list(
    label       = label,
    phi_scale   = phi_scale,
    k1_scale    = k1_scale,
    sensitivity = sens,
    scalars     = tibble::tibble(k1 = dd$k1, mean_log_phi = dd$mean_log_phi),
    log_phi_med = matrixStats::colMedians(lp),
    log_phi_sd  = matrixStats::colSds(lp),
    log_phi_draws = extract_log_phi_draws(
      draws_obj = dd, n_individuals = stan_data$n_individuals,
      thin_to = thin_to, seed = seed
    ),
    diag = fit$summary(c("k1", "mean_log_phi"), "rhat", "ess_bulk")
  )
}

ps_collect <- function(...) {
  fits <- list(...)
  stats::setNames(fits, vapply(fits, function(f) f$label, character(1)))[PS_LEVELS]
}

ps_shift_summary <- function(fits, virus, ref = "baseline") {
  base_med <- fits[[ref]]$log_phi_med
  between  <- stats::sd(base_med)
  
  purrr::map_dfr(fits, function(f) {
    dif   <- f$log_phi_med - base_med
    delta <- stats::median(dif)
    tibble::tibble(
      virus          = virus,
      label          = f$label,
      phi_scale      = f$phi_scale,
      k1_scale       = f$k1_scale,
      median_log_phi = stats::median(f$log_phi_med),
      delta          = delta,
      resid_sd       = stats::sd(dif - delta),
      resid_frac     = stats::sd(dif - delta) / between,
      spearman       = stats::cor(f$log_phi_med, base_med, method = "spearman"),
      k1_med         = stats::median(f$scalars$k1),
      max_rhat       = max(f$diag$rhat),
      min_ess        = min(f$diag$ess_bulk)
    )
  })
}

ps_sensitivity_table <- function(fits, virus) {
  purrr::map_dfr(fits, "sensitivity") |> dplyr::mutate(virus = virus, .before = 1)
}

# stage 2

fit_ps_stage2 <- function(model, sero_stan_data, label,
                          seed = 2, iter_warmup = 200, iter_sampling = 1000,
                          chains = 4, parallel_chains = 4, refresh = 200,
                          threads_per_chain = 1) {
  
  fit <- model$sample(
    data = sero_stan_data, seed = seed,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    chains = chains, parallel_chains = parallel_chains, refresh = refresh,
    threads_per_chain = threads_per_chain
  )
  
  dd      <- posterior::as_draws_df(fit$draws())
  age_max <- sero_stan_data$age_max
  
  lam <- as.matrix(dd[, paste0("lambda_long[", seq_len(age_max), "]"), drop = FALSE])
  pb  <- as.matrix(dd[, grep("^prob_by_group\\[", names(dd), value = TRUE), drop = FALSE])
  
  S <- t(apply(lam, 1, function(x) exp(-cumsum(x)))) # seroprevalence (w/out seroreversion)
  
  scalars <- tibble::tibble(
    lambda_1  = sr_get_draw(dd, "lambda_1"),
    kappa     = sr_get_draw(dd, "kappa"),
    psi       = sr_get_draw(dd, "psi"),
    sigma_min = matrixStats::rowMins(
      as.matrix(dd[, grep("^sigma_phi\\[", names(dd), value = TRUE), drop = FALSE]))
  )
  
  list(
    label   = label,
    age_max = age_max,
    scalars = scalars,
    foi = tibble::tibble(
      label  = label,
      Age    = seq_len(age_max),
      median = matrixStats::colMedians(lam),
      lo     = matrixStats::colQuantiles(lam, probs = 0.05),
      hi     = matrixStats::colQuantiles(lam, probs = 0.95)
    ),
    seroprev = tibble::tibble(
      label  = label,
      Age    = seq_len(age_max),
      median = 1 - matrixStats::colMedians(S)
    ),
    prob_mean = colMeans(pb),
    diag = fit$summary(c("lambda_1", "kappa", "sigma_phi"), "rhat", "ess_bulk")
  )
}

ps_collect_s2 <- function(...) {
  fits <- list(...)
  stats::setNames(fits, vapply(fits, function(f) f$label, character(1)))[PS_LEVELS]
}

# summaries
ps_foi_summary <- function(s2fits, virus, ref = "baseline",
                           ages = c(1, 2, 5, 10, 20, 40)) {
  base  <- s2fits[[ref]]
  pars  <- setdiff(names(base$scalars), "sigma_min")
  b_med <- vapply(base$scalars[pars], stats::median, numeric(1))
  b_sd  <- vapply(base$scalars[pars], stats::sd,     numeric(1))
  ages  <- ages[ages <= base$age_max]
  
  purrr::map_dfr(s2fits, function(f) {
    med <- vapply(f$scalars[pars], stats::median, numeric(1))
    sp  <- f$seroprev$median[match(ages, f$seroprev$Age)]
    bsp <- base$seroprev$median[match(ages, base$seroprev$Age)]
    dplyr::bind_rows(
      tibble::tibble(quantity = pars, value = med,
                     z = (med - b_med) / b_sd, abs_diff = med - b_med),
      tibble::tibble(quantity = paste0("seroprev_age", ages), value = sp,
                     z = NA_real_, abs_diff = sp - bsp)
    ) |>
      dplyr::mutate(virus = virus, label = f$label,
                    sigma_min = stats::median(f$scalars$sigma_min), .before = 1)
  })
}

ps_state_prob_summary <- function(s2fits, virus, ref = "baseline") {
  base <- s2fits[[ref]]$prob_mean
  purrr::map_dfr(s2fits, function(f) {
    tibble::tibble(
      virus         = virus,
      label         = f$label,
      max_abs_diff  = max(abs(f$prob_mean - base)),
      mean_abs_diff = mean(abs(f$prob_mean - base)),
      prop_gt_0.05  = mean(abs(f$prob_mean - base) > 0.05)
    )
  })
}

ps_reclass_summary <- function(s2fits, virus, k, ref = "baseline") {
  n     <- length(s2fits[[ref]]$prob_mean) / k
  stopifnot(n == round(n))
  modal <- function(p) apply(matrix(p, nrow = n, ncol = k), 1, which.max)
  base  <- modal(s2fits[[ref]]$prob_mean)
  
  purrr::map_dfr(s2fits, function(f) {
    m <- modal(f$prob_mean)
    tibble::tibble(virus = virus, label = f$label, n = n,
                   n_reclassified    = sum(m != base),
                   prop_reclassified = mean(m != base))
  })
}

# figure

ps_scale <- function() ggplot2::scale_colour_viridis_d(end = 0.9, option = "C")

plot_ps_logphi <- function(fits, shift_tab) {
  base_med <- fits[["baseline"]]$log_phi_med
  
  d <- purrr::map_dfr(fits, function(f) {
    tibble::tibble(label = f$label, base = base_med, rung = f$log_phi_med)
  }) |>
    dplyr::left_join(dplyr::select(shift_tab, label, delta), by = "label") |>
    dplyr::mutate(label = ps_label_factor(label)) |>
    dplyr::filter(label != "baseline") 
  
  ggplot2::ggplot(d, ggplot2::aes(base, rung, colour = label)) +
    ggplot2::geom_abline(ggplot2::aes(intercept = delta, slope = 1),
                         linetype = "dashed", linewidth = 0.5, colour = "grey30") +
    ggplot2::geom_point(alpha = 0.3, size = 0.5, show.legend = FALSE) +
    ggplot2::facet_wrap(~label, nrow = 1) +
    ps_scale() +
    ggplot2::labs(x = expression(log~phi[i]~"(baseline prior)"),
                  y = expression(log~phi[i]~"(varied prior)")) +
    ggplot2::theme_bw(base_size = 9)
}

plot_ps_state_probs <- function(s2fits) {
  base <- s2fits[["baseline"]]$prob_mean
  
  d <- purrr::map_dfr(s2fits, function(f) {
    tibble::tibble(label = f$label, base = base, rung = f$prob_mean)
  }) |>
    dplyr::mutate(label = ps_label_factor(label)) |>
    dplyr::filter(label != "baseline") 
  
  ggplot2::ggplot(d, ggplot2::aes(base, rung, colour = label)) +
    ggplot2::geom_abline(linetype = "dashed", linewidth = 0.5, colour = "grey30") +
    ggplot2::geom_point(alpha = 0.3, size = 0.5, show.legend = FALSE) +
    ggplot2::facet_wrap(~label, nrow = 1) +
    ps_scale() +
    ggplot2::labs(x = "P(serostate), baseline prior",
                  y = "P(serostate), varied prior") +
    ggplot2::theme_bw(base_size = 9)
}

plot_ps_foi <- function(s2fits) {
  d <- purrr::map_dfr(s2fits, "foi") |>
    dplyr::mutate(label = ps_label_factor(label))
  
  ggplot2::ggplot(d, ggplot2::aes(Age, median, colour = label, fill = label)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(linetype = label), linewidth = 0.65) +
    ps_scale() +
    ggplot2::scale_fill_viridis_d(end = 0.9, option = "C") +
    ggplot2::scale_linetype_manual(
      values = c("solid", "twodash", "longdash", "dashed", "dotdash", "dotted")) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Age (Years)", y = "Force of Infection",
                  colour = NULL, fill = NULL, linetype = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")
}

make_ps_figure <- function(fits, s2fits, shift_tab) {
  (plot_ps_logphi(fits, shift_tab) /
     plot_ps_state_probs(s2fits) /
     plot_ps_foi(s2fits)) +
    patchwork::plot_layout(heights = c(1, 1, 1.2)) +
    patchwork::plot_annotation(tag_levels = "A")
}

# supplementary table
make_ps_table <- function(shift_tab, foi_tab, state_tab, cens_tab, reclass_tab) {
  foi_wide <- foi_tab |>
    dplyr::filter(quantity %in% c("lambda_1", "kappa")) |>
    dplyr::select(virus, label, quantity, z) |>
    tidyr::pivot_wider(names_from = quantity, values_from = z,
                       names_prefix = "z_")
  
  tab <- shift_tab |>
    dplyr::left_join(foi_wide, by = c("virus", "label")) |>
    dplyr::left_join(dplyr::select(state_tab, virus, label, max_abs_diff),
                     by = c("virus", "label")) |>
    dplyr::left_join(dplyr::select(cens_tab, virus, prop_censored), by = "virus") |>
    dplyr::left_join(dplyr::select(reclass_tab, virus, label, prop_reclassified),
                     by = c("virus", "label")) |>
    dplyr::mutate(label = ps_label_factor(label)) |>
    dplyr::arrange(virus, label) |>
    dplyr::select(virus, label, phi_scale, k1_scale, delta, resid_sd, resid_frac,
                  spearman, k1_med, z_lambda_1, z_kappa, max_abs_diff,
                  prop_reclassified, prop_censored, max_rhat, min_ess)
  
  gt::gt(tab, groupname_col = "virus") |>
    gt::tab_header(
      title = "Sensitivity of Stage 1 and Stage 2 inference to the priors on phi and k1",
      subtitle = paste(
        "delta: common shift in median log phi relative to baseline.",
        "resid_frac: residual SD as a fraction of the between-individual SD.",
        "z: shift in FOI parameter median in baseline posterior SDs.")) |>
    gt::fmt_number(columns = c(delta, resid_sd, resid_frac, spearman, k1_med,
                               z_lambda_1, z_kappa, max_abs_diff, max_rhat),
                   decimals = 3) |>
    gt::fmt_number(columns = c(prop_censored, prop_reclassified), decimals = 3) |>
    gt::fmt_number(columns = min_ess, decimals = 0) |>
    gt::cols_label(label = "PRIOR", phi_scale = "PHI SCALE", k1_scale = "K1 SCALE",
                   delta = "DELTA", resid_sd = "RESID SD", resid_frac = "RESID FRAC",
                   spearman = "SPEARMAN", k1_med = "K1", z_lambda_1 = "Z LAMBDA1",
                   z_kappa = "Z KAPPA", max_abs_diff = "MAX |dP|", prop_reclassified = "PROP RECLASS",
                   prop_censored = "PROP CENS", max_rhat = "MAX RHAT",
                   min_ess = "MIN ESS")
}