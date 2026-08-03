#### adapt existing sero_stan_data list for revised SIS/SIIS Stan files ####

add_sero_switches <- function(stan_data,
                              estimate_omega = 1L,
                              omega_fixed    = 0,
                              omega_prior_sd = 0.1,
                              sigma_by_state = 0L,
                              foi_piecewise  = 0L,
                              bracket_cuts   = c(0, 5, 10, 20, 40)) {
  
  need <- c("n", "age_max", "ages", "n_phi_draws", "log_phi_draws")
  stopifnot(all(need %in% names(stan_data)))
  
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
    bracket_of_age = as.array(bracket_of_age)
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
  
  ## A named log_lik_var pins each model; an unnamed one is a candidate set.
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
  
  ## --- validity guard: loo_compare is meaningless across different data ------
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
## refits antibody_mech_k0.stan with k0 held at -2, -1, 0, 1, 2 [and with k0 freely estimated]
##
## panel  image
##
##   A    fitted dose-response curves
##   B    log(phi)
##   C    FOI posterior
##

K0_SPECS <- c("-2", "-1", "0", "1", "2", "free")
K0_IDS   <- c("m2", "m1", "z0", "p1", "p2", "free")
K0_FIXED <- c("-2", "-1", "0", "1", "2")
K0_PHI_PRIOR_SCALE <- 500
K0_PRIOR_SD        <- 5

K0_S2_THIN     <- 500 #NULL #100
K0_S2_CHAINS   <- 4
K0_S2_WARMUP   <- 200
K0_S2_SAMPLING <- 800


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
  cmdstanr::cmdstan_model(stan_file, dir = dir)
}

fit_k0_stage1 <- function(model, stan_data, spec,
                          keep_log_phi_draws = FALSE, thin_to = NULL,
                          seed = 1, iter_warmup = 500, iter_sampling = 500,
                          chains = 4, parallel_chains = 4, refresh = 200) {
  
  free <- identical(spec, "free")
  zeta <- if (free) 0 else as.numeric(spec)
  
  stan_data$estimate_k0     <- as.integer(free)
  stan_data$transform_prior <- as.integer(!free)
  stan_data$k0_fixed        <- zeta
  stan_data$k0_prior_sd     <- K0_PRIOR_SD
  stan_data$phi_prior_scale <- K0_PHI_PRIOR_SCALE
  
  fit <- model$sample(
    data = stan_data, seed = seed,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    chains = chains, parallel_chains = parallel_chains, refresh = refresh
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
                          dir = "Stan/bin", seed = 2,
                          iter_warmup = K0_S2_WARMUP,
                          iter_sampling = K0_S2_SAMPLING,
                          chains = K0_S2_CHAINS, parallel_chains = K0_S2_CHAINS,
                          refresh = 0) {
  
  ensure_dir(dir)
  model <- cmdstanr::cmdstan_model(sero_stan_file, dir = dir)
  
  fit <- model$sample(
    data = sero_stan_data, seed = seed,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    chains = chains, parallel_chains = parallel_chains, refresh = refresh
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

## mixture median of phi at a given age
##
## solve sum_x w_x Phi((log m - mu_x)/sigma) = 0.5

mixture_phi_median <- function(w, mu_x, sigma) {
  w <- w[seq_along(mu_x)]
  if (sum(w) <= 0) return(NA_real_)
  w <- w / sum(w)
  f <- function(z) sum(w * stats::pnorm((z - mu_x) / sigma)) - 0.5
  exp(stats::uniroot(f, lower = min(mu_x) - 6 * sigma,
                     upper = max(mu_x) + 6 * sigma)$root)
}

#' Mean and median of the fitted mixture at every age
mixture_phi_summaries <- function(probs_age, mu_x, sigma) {
  k <- length(mu_x)
  P <- probs_age[, seq_len(k), drop = FALSE]
  data.frame(
    Age    = seq_len(nrow(P)) - 1L,
    mean   = as.numeric(P %*% exp(mu_x + 0.5 * sigma^2)),
    median = apply(P, 1, mixture_phi_median, mu_x = mu_x, sigma = sigma)
  )
}