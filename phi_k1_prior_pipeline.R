## prior sensitivity analysis on phi and k1 

PS_VIRUSES <- c("CVA6", "EV71", "EV68")

PS_SERO_FILE <- c(
  CVA6 = "Stan/antibody_sero.stan",                  # SII
  EV71 = "Stan/antibody_sero.stan",                  # SII
  EV68 = "Stan/antibody_sero_two_compartment.stan"   # SI
)

# tight/wide priors for phi and k1 + replicate for MCMC noise floor
PS_SPECS <- tibble::tribble(
  ~label,      ~id,     ~phi_scale, ~k1_scale, ~seed_offset,
  "baseline",  "base",         500,      1.00,           0L,
  "replicate", "rep",          500,      1.00,           1L,
  "phi tight", "phiT",          50,      1.00,           0L,
  "phi wide",  "phiW",        5000,      1.00,           0L,
  "k1 tight",  "k1T",          500,      0.25,           0L,
  "k1 wide",   "k1W",          500,      4.00,           0L
)

PS_N_STATES <- c(CVA6 = 3L, EV71 = 3L, EV68 = 2L)

## grids

ps_grid <- tidyr::expand_grid(virus = PS_VIRUSES, PS_SPECS) |>
  dplyr::mutate(prep = rlang::syms(paste0("ab_prep_", virus)))

ps_grid_s2 <- ps_grid |>
  dplyr::mutate(
    fit        = rlang::syms(paste0("ps_fit_", virus, "_", id)),
    sero_model = rlang::syms(ifelse(virus == "EV68",
                                    "ps_sero_model_si", "ps_sero_model_sii"))
  )

## targets

ps_targets <- list(

  tar_target(ps_sero_file_sii, "Stan/antibody_sero.stan", format = "file"),
  tar_target(ps_sero_file_si,  "Stan/antibody_sero_two_compartment.stan",
             format = "file"),
  tar_target(ps_sero_model_sii,
             cmdstanr::cmdstan_model(ps_sero_file_sii, dir = "Stan/bin",
                                     cpp_options = list(stan_threads = TRUE))),
  tar_target(ps_sero_model_si,
             cmdstanr::cmdstan_model(ps_sero_file_si, dir = "Stan/bin",
                                     cpp_options = list(stan_threads = TRUE))),

  tar_target(ps_cens_table, dplyr::bind_rows(
    ps_censoring_summary(ab_prep_CVA6, "CVA6"),
    ps_censoring_summary(ab_prep_EV71, "EV71"),
    ps_censoring_summary(ab_prep_EV68, "EV68"))),

  ## Stage 1
  tar_map(
    values = ps_grid,
    names = c(virus, id),
    tar_target(
      ps_fit,
      fit_ps_stage1(
        model     = k0_model,            # same compiled binary as k0 analysis
        stan_data = prep$stan_data,
        label     = label,
        phi_scale = phi_scale,
        k1_scale  = k1_scale,
        thin_to   = n_phi_draws_thin,
        seed              = ab_seed + seed_offset,
        iter_warmup       = ab_iter_warmup,
        iter_sampling     = ab_iter_sampling,
        chains            = ab_chains,
        parallel_chains   = ab_parallel,
        refresh           = ab_refresh,
        threads_per_chain = ab_threads
      )
    )
  ),

  tar_target(ps_fits_CVA6, ps_collect(ps_fit_CVA6_base, ps_fit_CVA6_rep,
                                      ps_fit_CVA6_phiT, ps_fit_CVA6_phiW, 
                                      ps_fit_CVA6_k1T,  ps_fit_CVA6_k1W)),
  tar_target(ps_fits_EV71, ps_collect(ps_fit_EV71_base, ps_fit_EV71_rep,
                                      ps_fit_EV71_phiT, ps_fit_EV71_phiW, 
                                      ps_fit_EV71_k1T,  ps_fit_EV71_k1W)),
  tar_target(ps_fits_EV68, ps_collect(ps_fit_EV68_base, ps_fit_EV68_rep,
                                      ps_fit_EV68_phiT,  ps_fit_EV68_phiW, 
                                      ps_fit_EV68_k1T,  ps_fit_EV68_k1W)),

  ## Stage 2 
  tar_map(
    values = ps_grid_s2,
    names = c(virus, id),
    tar_target(
      ps_sero_data,
      make_sero_stan_data(
        serum_meta     = prep$serum_meta,
        log_phi_draws  = fit$log_phi_draws,
        sigma_by_state = 0
      )
    ),
    tar_target(
      ps_foi,
      fit_ps_stage2(
        model          = sero_model,
        sero_stan_data = ps_sero_data,
        label          = label,
        seed              = sero_seed + seed_offset,
        iter_warmup       = sero_iter_warmup,
        iter_sampling     = sero_iter_sampling - sero_iter_warmup,
        chains            = sero_chains,
        parallel_chains   = sero_parallel,
        refresh           = sero_refresh,
        threads_per_chain = sero_threads
      )
    )
  ),

  tar_target(ps_s2_CVA6, ps_collect_s2(ps_foi_CVA6_base, ps_foi_CVA6_rep,
                                       ps_foi_CVA6_phiT, ps_foi_CVA6_phiW, 
                                       ps_foi_CVA6_k1T,  ps_foi_CVA6_k1W)),
  tar_target(ps_s2_EV71, ps_collect_s2(ps_foi_EV71_base, ps_foi_EV71_rep,
                                       ps_foi_EV71_phiT, ps_foi_EV71_phiW, 
                                       ps_foi_EV71_k1T,  ps_foi_EV71_k1W)),
  tar_target(ps_s2_EV68, ps_collect_s2(ps_foi_EV68_base, ps_foi_EV68_rep,
                                       ps_foi_EV68_phiT, ps_foi_EV68_phiW, 
                                       ps_foi_EV68_k1T,  ps_foi_EV68_k1W)),

  ## summary tables
  tar_target(ps_shift_table, dplyr::bind_rows(
    ps_shift_summary(ps_fits_CVA6, "CVA6"),
    ps_shift_summary(ps_fits_EV71, "EV71"),
    ps_shift_summary(ps_fits_EV68, "EV68"))),

  tar_target(ps_sens_table, dplyr::bind_rows(
    ps_sensitivity_table(ps_fits_CVA6, "CVA6"),
    ps_sensitivity_table(ps_fits_EV71, "EV71"),
    ps_sensitivity_table(ps_fits_EV68, "EV68"))),

  tar_target(ps_foi_table, dplyr::bind_rows(
    ps_foi_summary(ps_s2_CVA6, "CVA6"),
    ps_foi_summary(ps_s2_EV71, "EV71"),
    ps_foi_summary(ps_s2_EV68, "EV68"))),

  tar_target(ps_state_table, dplyr::bind_rows(
    ps_state_prob_summary(ps_s2_CVA6, "CVA6"),
    ps_state_prob_summary(ps_s2_EV71, "EV71"),
    ps_state_prob_summary(ps_s2_EV68, "EV68"))),
  
  tar_target(ps_reclass_table, dplyr::bind_rows(
    ps_reclass_summary(ps_s2_CVA6, "CVA6", PS_N_STATES[["CVA6"]]),
    ps_reclass_summary(ps_s2_EV71, "EV71", PS_N_STATES[["EV71"]]),
    ps_reclass_summary(ps_s2_EV68, "EV68", PS_N_STATES[["EV68"]]))),

  ## figures
  tar_target(
    file_ps_CVA6,
    save_plot_pdf(
      make_ps_figure(ps_fits_CVA6, ps_s2_CVA6,
                     dplyr::filter(ps_shift_table, virus == "CVA6")),
      "outputs/CVA6/CVA6_prior_sensitivity.pdf", 11, 9),
    format = "file"
  ),
  tar_target(
    file_ps_EV71,
    save_plot_pdf(
      make_ps_figure(ps_fits_EV71, ps_s2_EV71,
                     dplyr::filter(ps_shift_table, virus == "EV71")),
      "outputs/EV71/EV71_prior_sensitivity.pdf", 11, 9),
    format = "file"
  ),
  tar_target(
    file_ps_EV68,
    save_plot_pdf(
      make_ps_figure(ps_fits_EV68, ps_s2_EV68,
                     dplyr::filter(ps_shift_table, virus == "EV68")),
      "outputs/EV68/EV68_prior_sensitivity.pdf", 11, 9),
    format = "file"
  ),

  ## supplementary table
  tar_target(
    file_ps_table,
    save_gt_pdf(
      make_ps_table(ps_shift_table, ps_foi_table, ps_state_table,
                    ps_cens_table, ps_reclass_table),
      "outputs/prior_sensitivity/prior_sensitivity_table.pdf"),
    format = "file"
  )
)