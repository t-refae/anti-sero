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
                             omega_prior_sd = 0.1) {
  
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