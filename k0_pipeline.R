## k0 sensitivity analysis

K0_VIRUSES <- c("CVA6")

k0_grid <- tidyr::expand_grid(
  virus = K0_VIRUSES,
  spec  = K0_SPECS
) |>
  dplyr::mutate(
    id   = K0_IDS[match(spec, K0_SPECS)],
    prep = rlang::syms(paste0("ab_prep_", virus)),
    keep_draws = virus == "CVA6" & spec %in% K0_FIXED
  )

k0_grid_s2 <- k0_grid |>
  dplyr::filter(keep_draws) |>
  dplyr::mutate(fit = rlang::syms(paste0("k0_fit_", virus, "_", id)))

k0_targets <- list(
  
  tar_target(k0_model_file, "Stan/antibody_mech_k0.stan", format = "file"),
  tar_target(k0_model, k0_compile(k0_model_file)),
  
  ## antibody fitting
  tar_map(
    values = k0_grid,
    names = c(virus, id),
    tar_target(
      k0_fit,
      fit_k0_stage1(
        model     = k0_model,
        stan_data = prep$stan_data,
        spec      = spec,
        keep_log_phi_draws = keep_draws,
        thin_to   = K0_S2_THIN,
        seed      = ab_seed,
        iter_warmup   = ab_iter_warmup,
        iter_sampling = ab_iter_sampling,
        chains        = ab_chains,
        parallel_chains = ab_parallel,
        refresh   = ab_refresh
      )
    )
  ),
  
  ## serocatalytic fitting
  tar_map(
    values = k0_grid_s2,
    names = c(virus, id),
    tar_target(
      k0_sero_data,
      make_sero_stan_data(
        serum_meta    = prep$serum_meta,
        log_phi_draws = fit$log_phi_draws
      )
    ),
    tar_target(
      k0_foi,
      fit_k0_stage2(
        sero_stan_file = "Stan/antibody_sero.stan",   # SII, as used for CVA6
        sero_stan_data = k0_sero_data,
        spec           = spec,
        seed           = sero_seed,
        # stan_cpp=stan_cpp,
        threads_per_chain=sero_threads
      )
    )
  ),
  
  tar_target(
    k0_fits_CVA6,
    k0_collect(k0_fit_CVA6_m2, k0_fit_CVA6_m1, k0_fit_CVA6_z0,
               k0_fit_CVA6_p1, k0_fit_CVA6_p2, k0_fit_CVA6_free)
  ),
  tar_target(
    k0_foi_CVA6,
    dplyr::bind_rows(k0_foi_CVA6_m2, k0_foi_CVA6_m1, k0_foi_CVA6_z0,
                     k0_foi_CVA6_p1, k0_foi_CVA6_p2)
  ),
  
  tar_target(k0_plot_CVA6, make_k0_figure_cva6(k0_fits_CVA6, k0_foi_CVA6)),
  tar_target(
    file_k0_CVA6,
    save_plot_pdf(k0_plot_CVA6, "outputs/CVA6/CVA6_k0_invariance.pdf", 9, 7),
    format = "file"
  ),
  
  ## Numbers for the caption and the supplementary table ----------------------
  tar_target(k0_curve_diff_CVA6, k0_max_curve_diff(k0_fits_CVA6[K0_FIXED])),
  tar_target(k0_table_df, make_k0_table(k0_fits_CVA6, "CVA6")),
  tar_target(
    file_k0_table_csv,
    {
      ensure_dir("outputs/CVA6")
      write.csv(k0_table_df, "outputs/CVA6/CVA6_k0_invariance.csv",
                row.names = FALSE)
      "outputs/CVA6/CVA6_k0_invariance.csv"
    },
    format = "file"
  )
)

