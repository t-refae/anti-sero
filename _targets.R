# _targets.R
library(targets)
library(stantargets)
library(tarchetypes)

# global options
tar_option_set(
  packages = c(
    # core functionality
    "targets", "stantargets", "tarchetypes", "cmdstanr",  "loo",
    
    # data analysis tools
    "dplyr", "tidyr", "purrr", "stringr", "posterior", "matrixStats", "reshape2",
    
    # plotting / saving pdfs
    "ggplot2", "gt", "gtExtras", "cowplot", "grid", "gridExtra", "khroma",
    "patchwork", "scales", "pdftools", "showtext", "viridis", "GGally"
  ),
  format = "qs"
)

# load all functions in R/ directory
tar_source()
source("k0_pipeline.R")

# to render tables
ensure_chromote_browser()

# ---------- Global MCMC settings ----------

# antibody model
ab_iter_warmup   <- 500
ab_iter_sampling <- 500
ab_chains        <- 4
ab_parallel      <- 4
ab_seed          <- 1
ab_refresh       <- 200

# serocatalytic model
sero_iter_sampling <- 1000
sero_iter_warmup   <- floor(sero_iter_sampling / 5)
sero_chains        <- 4
sero_parallel      <- 4
sero_seed          <- 2
sero_refresh       <- 200

# seroreversion fits: omega is weakly identified and needs longer adaptation
sero_sr_iter_warmup <- 800
sero_sr_adapt_delta <- 0.9
sero_omega_prior_sd <- 0.5
sero_sigma_by_state <- 1L

# thin antibody posterior draws for speed (NULL = no thinning)
n_phi_draws_thin <- 500 #NULL 

# ---------- Pipeline ----------
list(
  # -------- Antibody data prep (raw well observations -> stan_data + serum_meta) --------
  tar_target(ab_prep_CVA6, prep_antibody_data(virus = "CVA6", raw_prefix = "CA6")),
  tar_target(ab_prep_EV71, prep_antibody_data(virus = "EV71", raw_prefix = "EV71")),
  tar_target(ab_prep_EV68, prep_antibody_data(virus = "EV68", raw_prefix = "E68")),
  
  # -------- Antibody fits (mechanistic model) --------
  tar_stan_mcmc(
    name = "ab_CVA6",
    stan_files = "Stan/antibody_mech.stan",
    data = ab_prep_CVA6$stan_data,
    iter_warmup = ab_iter_warmup,
    iter_sampling = ab_iter_sampling,
    chains = ab_chains,
    parallel_chains = ab_parallel,
    seed = ab_seed,
    refresh = ab_refresh
  ),
  
  tar_stan_mcmc(
    name = "ab_EV71",
    stan_files = "Stan/antibody_mech.stan",
    data = ab_prep_EV71$stan_data,
    iter_warmup = ab_iter_warmup,
    iter_sampling = ab_iter_sampling,
    chains = ab_chains,
    parallel_chains = ab_parallel,
    seed = ab_seed,
    refresh = ab_refresh
  ),
  
  tar_stan_mcmc(
    name = "ab_EV68",
    stan_files = "Stan/antibody_mech.stan",
    data = ab_prep_EV68$stan_data,
    iter_warmup = ab_iter_warmup,
    iter_sampling = ab_iter_sampling,
    chains = ab_chains,
    parallel_chains = ab_parallel,
    seed = ab_seed,
    refresh = ab_refresh
  ),
  
  # --- Antibody raw data for PPC  ---
  tar_target(ab_raw_CVA6, read.csv("data/raw/CA6_raw_well_observations.csv")),
  tar_target(ab_raw_EV71, read.csv("data/raw/EV71_raw_well_observations.csv")),
  tar_target(ab_raw_EV68, read.csv("data/raw/E68_raw_well_observations.csv")),
  
  # --- Prior + Posterior PPC overlay plots ---
  tar_target(
    ab_ppc_CVA6,
    antibody_ppc_overlay_plot(
      virus_prefix = "CVA6",
      raw_df = ab_raw_CVA6,
      draws_obj = ab_CVA6_draws_antibody_mech,
      thin_draws_post = 1000,
      ncol = 4,
      n_prior_lines = 500,
      prior_alpha = 0.08,      
      prior_linewidth = 0.25
    )
  ),
  
  tar_target(
    ab_ppc_EV71,
    antibody_ppc_overlay_plot(
      virus_prefix = "EV71",
      raw_df = ab_raw_EV71,
      draws_obj = ab_EV71_draws_antibody_mech,
      thin_draws_post = 1000,
      ncol = 4,
      n_prior_lines = 500,        
      prior_alpha = 0.08,         
      prior_linewidth = 0.25
    )
  ),
  
  tar_target(
    ab_ppc_EV68,
    antibody_ppc_overlay_plot(
      virus_prefix = "EV68",
      raw_df = ab_raw_EV68,
      draws_obj = ab_EV68_draws_antibody_mech,
      thin_draws_post = 1000,
      ncol = 4,
      n_prior_lines = 500,
      prior_alpha = 0.08,
      prior_linewidth = 0.25
    )
  ),
  
  # -------- Extract posterior draws of log(phi) per individual (with optional thinning) --------
  tar_target(
    log_phi_draws_CVA6,
    extract_log_phi_draws(
      draws_obj = ab_CVA6_draws_antibody_mech,
      n_individuals = ab_prep_CVA6$stan_data$n_individuals,
      thin_to = n_phi_draws_thin,
      seed = 1
    )
  ),
  
  tar_target(
    log_phi_draws_EV71,
    extract_log_phi_draws(
      draws_obj = ab_EV71_draws_antibody_mech,
      n_individuals = ab_prep_EV71$stan_data$n_individuals,
      thin_to = n_phi_draws_thin,
      seed = 1
    )
  ),
  
  tar_target(
    log_phi_draws_EV68,
    extract_log_phi_draws(
      draws_obj = ab_EV68_draws_antibody_mech,
      n_individuals = ab_prep_EV68$stan_data$n_individuals,
      thin_to = n_phi_draws_thin,
      seed = 1
    )
  ),
  
  # -------- Build serocatalytic Stan data --------
  tar_target(
    sero_stan_data_CVA6,
    make_sero_stan_data(
      serum_meta = ab_prep_CVA6$serum_meta,
      log_phi_draws = log_phi_draws_CVA6,
      sigma_by_state = sero_sigma_by_state
    )
  ),
  
  tar_target(
    sero_stan_data_EV71,
    make_sero_stan_data(
      serum_meta = ab_prep_EV71$serum_meta,
      log_phi_draws = log_phi_draws_EV71,
      sigma_by_state = sero_sigma_by_state
    )
  ),
  
  tar_target(
    sero_stan_data_EV68,
    make_sero_stan_data(
      serum_meta = ab_prep_EV68$serum_meta,
      log_phi_draws = log_phi_draws_EV68,
      sigma_by_state = sero_sigma_by_state
    )
  ),
  
  # -------- Build plotting data --------
  tar_target(
    sero_plot_data_CVA6,
    make_sero_plot_data(
      serum_meta = ab_prep_CVA6$serum_meta,
      log_phi_draws = log_phi_draws_CVA6
    )
  ),
  
  tar_target(
    sero_plot_data_EV71,
    make_sero_plot_data(
      serum_meta = ab_prep_EV71$serum_meta,
      log_phi_draws = log_phi_draws_EV71
    )
  ),

  tar_target(
    sero_plot_data_EV68,
    make_sero_plot_data(
      serum_meta = ab_prep_EV68$serum_meta,
      log_phi_draws = log_phi_draws_EV68
    )
  ),
  
  # -------- Serocatalytic fits --------
  tar_stan_mcmc( # SII model
    name = "CVA6",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_stan_data_CVA6,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  tar_stan_mcmc( # SII model
    name = "EV71",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_stan_data_EV71,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  tar_stan_mcmc( # SI model
    name = "EV68",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_stan_data_EV68,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  tar_stan_mcmc( # SII model
    name = "EV68_three_comp",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_stan_data_EV68,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  # -------- Combine draws & summaries --------
  tar_target(
    draws_combined,
    combine_draws_three_viruses(
      CVA6_draws_antibody_sero,
      EV71_draws_antibody_sero,
      EV68_draws_antibody_sero_two_compartment
    )
  ),
  
  tar_target(
    summaries_combined,
    combine_summaries_three_viruses(
      CVA6_summary_antibody_sero,
      EV71_summary_antibody_sero,
      EV68_summary_antibody_sero_two_compartment
    )
  ),
  
  # -------- Phi plots --------
  tar_target(
    CVA6_phi_plot,
    plot_phi_points("CVA6", sero_plot_data_CVA6, CVA6_summary_antibody_sero)
  ),
  tar_target(
    EV71_phi_plot,
    plot_phi_points("EV71", sero_plot_data_EV71, EV71_summary_antibody_sero)
  ),
  tar_target(
    EV68_phi_plot,
    plot_phi_points("EV68", sero_plot_data_EV68, EV68_summary_antibody_sero_two_compartment)
  ),
  
  # -------- Weighted-average phi plot --------
  tar_target(
    CVA6_weighted_avg_phi_plot,
    plot_weighted_avg_phi("CVA6", sero_plot_data_CVA6,
                          CVA6_draws_antibody_sero,
                          CVA6_summary_antibody_sero)
  ),
  tar_target(
    EV71_weighted_avg_phi_plot,
    plot_weighted_avg_phi("EV71", sero_plot_data_EV71,
                          EV71_draws_antibody_sero,
                          EV71_summary_antibody_sero)
  ),
  tar_target(
    EV68_weighted_avg_phi_plot,
    plot_weighted_avg_phi("EV68", sero_plot_data_EV68,
                          EV68_draws_antibody_sero_two_compartment,
                          EV68_summary_antibody_sero_two_compartment)
  ),
  
  # -------- Serodynamics plots --------
  tar_target(
    CVA6_serodynamics_plot,
    plot_serodynamics("CVA6", sero_plot_data_CVA6, CVA6_draws_antibody_sero)
  ),
  tar_target(
    EV71_serodynamics_plot,
    plot_serodynamics("EV71", sero_plot_data_EV71, EV71_draws_antibody_sero)
  ),
  tar_target( # SI model
    EV68_serodynamics_plot,
    plot_serodynamics("EV68", sero_plot_data_EV68, EV68_draws_antibody_sero_two_compartment)
  ),
  
  tar_target( # SII model
    EV68_serodynamics_plot_three_comp,
    plot_serodynamics("EV68", sero_plot_data_EV68, EV68_three_comp_draws_antibody_sero)
  ),
  
  # -------- FOI plots --------
  tar_target(
    CVA6_FOI_plot,
    plot_foi_short_long("CVA6", sero_plot_data_CVA6, CVA6_summary_antibody_sero,
                        CVA6_draws_antibody_sero)
  ),
  tar_target(
    EV71_FOI_plot,
    plot_foi_short_long("EV71", sero_plot_data_EV71, EV71_summary_antibody_sero,
                        EV71_draws_antibody_sero)
  ),
  tar_target(
    EV68_FOI_plot,
    plot_foi_short_long("EV68", sero_plot_data_EV68, EV68_summary_antibody_sero_two_compartment,
                        EV68_draws_antibody_sero_two_compartment)
  ),
  
  # single plot
  tar_target(
    combined_FOI_plot,
    plot_foi_combined(draws_combined, max_age = max(sero_plot_data_CVA6$age_max,
                                                    sero_plot_data_EV71$age_max,
                                                    sero_plot_data_EV68$age_max))
  ),
  
  # grid of plots
  tar_target(
    combined_FOI_grid_plot,
    plot_patchwork(
      CVA6_FOI_new_old_plot,
      EV71_FOI_new_old_plot,
      EV68_FOI_new_old_plot,
      FOI=TRUE
    )
  ),
  
  # new vs old FOI plots
  tar_target(
    CVA6_FOI_new_old_plot,
    plot_foi_new_old("CVA6", sero_plot_data_CVA6,
                     CVA6_summary_antibody_sero,
                     CVA6_draws_antibody_sero)
  ),
  tar_target(
    EV71_FOI_new_old_plot,
    plot_foi_new_old("EV71", sero_plot_data_EV71,
                     EV71_summary_antibody_sero,
                     EV71_draws_antibody_sero)
  ),
  tar_target(
    EV68_FOI_new_old_plot,
    plot_foi_new_old("EV68", sero_plot_data_EV68,
                     EV68_summary_antibody_sero_two_compartment,
                     EV68_draws_antibody_sero_two_compartment)
  ),
  
  
  # -------- Posterior density plots --------
  tar_target(
    combined_density_plots,
    plot_param_densities_combined(draws_combined)
  ),
  
  # -------- Parameter tables --------
  tar_target(
    CVA6_param_table,
    make_param_table("CVA6", CVA6_summary_antibody_sero)
  ),
  tar_target(
    EV71_param_table,
    make_param_table("EV71", EV71_summary_antibody_sero)
  ),
  tar_target(
    EV68_param_table,
    make_param_table("EV68", EV68_summary_antibody_sero_two_compartment)
  ),
  
  # -------- Save outputs --------
  tar_target(file_ab_ppc_CVA6, save_plot_pdf(ab_ppc_CVA6, "outputs/CVA6/CVA6_antibody_ppc_overlay.pdf", width = 12, height = 8), format = "file"),
  tar_target(file_ab_ppc_EV71, save_plot_pdf(ab_ppc_EV71, "outputs/EV71/EV71_antibody_ppc_overlay.pdf", width = 12, height = 8), format = "file"),
  tar_target(file_ab_ppc_EV68, save_plot_pdf(ab_ppc_EV68, "outputs/EV68/EV68_antibody_ppc_overlay.pdf", width = 12, height = 8), format = "file"),
  
  tar_target(file_CVA6_FOI, save_plot_pdf(CVA6_FOI_plot, "outputs/CVA6/CVA6_FOI_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV71_FOI, save_plot_pdf(EV71_FOI_plot, "outputs/EV71/EV71_FOI_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_FOI, save_plot_pdf(EV68_FOI_plot, "outputs/EV68/EV68_FOI_plot.pdf", 8, 6), format = "file"),
  
  tar_target(file_combined_FOI, save_plot_pdf(combined_FOI_plot, "outputs/combined/combined_FOI_plot.pdf", 8, 6), format = "file"),
  tar_target(file_combined_FOI_grid, save_plot_pdf(combined_FOI_grid_plot, "outputs/combined/combined_FOI_grid_plot.pdf", 12, 10), format = "file"),
  
  tar_target(file_density, save_plot_pdf(combined_density_plots, "outputs/combined/combined_density_plots.pdf", 10, 7), format = "file"),
  
  tar_target(file_CVA6_phi, save_plot_pdf(CVA6_phi_plot, "outputs/CVA6/CVA6_phi_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV71_phi, save_plot_pdf(EV71_phi_plot, "outputs/EV71/EV71_phi_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_phi, save_plot_pdf(EV68_phi_plot, "outputs/EV68/EV68_phi_plot.pdf", 8, 6), format = "file"),
  
  tar_target(file_CVA6_serodynamics, save_plot_pdf(CVA6_serodynamics_plot, "outputs/CVA6/CVA6_serodynamics_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV71_serodynamics, save_plot_pdf(EV71_serodynamics_plot, "outputs/EV71/EV71_serodynamics_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_serodynamics, save_plot_pdf(EV68_serodynamics_plot, "outputs/EV68/EV68_serodynamics_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_serodynamics_three_comp, save_plot_pdf(EV68_serodynamics_plot_three_comp, "outputs/EV68/EV68_serodynamics_plot_three_comp.pdf", 8, 6), format = "file"),
  
  tar_target(
    combined_serodynamics_plot,
    plot_patchwork(
      CVA6_serodynamics_plot,
      EV71_serodynamics_plot,
      EV68_serodynamics_plot
    )
  ),
  
  tar_target(file_combined_serodynamics, save_plot_pdf(combined_serodynamics_plot, "outputs/combined/combined_serodynamics_plot.pdf", 12, 10), format = "file"),
  
  tar_target(file_CVA6_weighted, save_plot_pdf(CVA6_weighted_avg_phi_plot, "outputs/CVA6/CVA6_weighted_avg_phi_plot.pdf", 10, 7), format = "file"),
  tar_target(file_EV71_weighted, save_plot_pdf(EV71_weighted_avg_phi_plot, "outputs/EV71/EV71_weighted_avg_phi_plot.pdf", 10, 7), format = "file"),
  tar_target(file_EV68_weighted, save_plot_pdf(EV68_weighted_avg_phi_plot, "outputs/EV68/EV68_weighted_avg_phi_plot.pdf", 10, 7), format = "file"),
  
  tar_target(
    combined_weighted_avg_phi_plot,
    plot_patchwork(
      CVA6_weighted_avg_phi_plot,
      EV71_weighted_avg_phi_plot,
      EV68_weighted_avg_phi_plot,
      lines=FALSE
    )
  ),
  
  tar_target(file_combined_weighted_avg_phi, save_plot_pdf(combined_weighted_avg_phi_plot, "outputs/combined/combined_weighted_avg_phi_plot.pdf", 12, 10), format = "file"),
  
  tar_target(file_CVA6_table, save_gt_pdf(CVA6_param_table, "outputs/CVA6/CVA6_param_table.pdf"), format = "file"),
  tar_target(file_EV71_table, save_gt_pdf(EV71_param_table, "outputs/EV71/EV71_param_table.pdf"), format = "file"),
  tar_target(file_EV68_table, save_gt_pdf(EV68_param_table, "outputs/EV68/EV68_param_table.pdf"), format = "file"),
  
  # save as concatenated pdfs (3 pages)
  tar_target(
    file_combined_parameter_tables,
    combine_pdfs(
      inputs = c(file_CVA6_table, file_EV71_table, file_EV68_table),
      output = "outputs/combined/combined_parameter_tables.pdf"
    ),
    format = "file"
  ),
  
  
  # -------- Reed Muench --------
  tar_target(
    reed_muench_data,
    make_reed_muench_data()
  ),
  
  tar_stan_mcmc(
    name = "ab_reed_muench",
    stan_files = "Stan/antibody_mech.stan",
    data = reed_muench_data$stan_data,
    iter_warmup = 1000,
    iter_sampling = 1000,
    chains = 4,
    parallel_chains = 4,
    seed = 1,
    refresh = 200
  ),
  
  tar_target(
    reed_muench_summary,
    summarise_reed_muench_fit(ab_reed_muench_draws_antibody_mech)
  ),
  
  tar_target(
    reed_muench_plot,
    plot_reed_muench_fit(
      raw_df = reed_muench_data$raw_df,
      draws_obj = ab_reed_muench_draws_antibody_mech
    )
  ),
  
  tar_target(file_reed_muench_plot, save_plot_pdf(reed_muench_plot, "outputs/reed_muench/reed_muench_fit.pdf", 8, 6), format = "file"),
  
  # -------- LOOCV --------
  tar_target(
    EV68_loo_compare,
    compare_si_sii_loo(
      si_draws = EV68_draws_antibody_sero_two_compartment,
      sii_draws = EV68_three_comp_draws_antibody_sero,
      virus = "EV68"
    )
  ),
  
  tar_target(
    file_EV68_loo_compare,
    {
      ensure_dir("outputs/EV68")
      write.csv(
        EV68_loo_compare,
        "outputs/EV68/EV68_SI_vs_SII_LOOCV.csv",
        row.names = FALSE
      )
      "outputs/EV68/EV68_SI_vs_SII_LOOCV.csv"
    },
    format = "file"
  ),
  
  tar_stan_mcmc(
    name = "CVA6_two_comp",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_stan_data_CVA6,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  tar_stan_mcmc(
    name = "EV71_two_comp",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_stan_data_EV71,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  tar_target(
    CVA6_loo_compare,
    compare_si_sii_loo(
      si_draws = CVA6_two_comp_draws_antibody_sero_two_compartment,
      sii_draws = CVA6_draws_antibody_sero,
      virus = "CVA6"
    )
  ),
  
  tar_target(
    EV71_loo_compare,
    compare_si_sii_loo(
      si_draws = EV71_two_comp_draws_antibody_sero_two_compartment,
      sii_draws = EV71_draws_antibody_sero,
      virus = "EV71"
    )
  ),
  
  tar_target(
    combined_loo_table_df,
    combine_loo_tables(
      cva6_loo = CVA6_loo_compare,
      ev71_loo = EV71_loo_compare,
      ev68_loo = EV68_loo_compare
    )
  ),
  
  tar_target(
    combined_loo_table,
    make_loo_table(combined_loo_table_df)
  ),
  
  tar_target(file_combined_loo_table, save_gt_pdf(combined_loo_table, "outputs/combined/combined_LOOCV_table.pdf"), format = "file"),
  
  tar_target(
    file_combined_loo_csv,
    {
      ensure_dir("outputs/combined")
      write.csv(
        combined_loo_table_df,
        "outputs/combined/combined_LOOCV_table.csv",
        row.names = FALSE
      )
      "outputs/combined/combined_LOOCV_table.csv"
    },
    format = "file"
  ),
  
  # -------- Misc. --------
  
  # SI Fig 1
  tar_target(
    survival_plots,
    plot_survival(
      baseline_k1 = 2,
      baseline_phi = 35
    )
  ),
  
  tar_target(file_survival_plots, save_plot_pdf(survival_plots, "outputs/misc/survival_plots.pdf", 8, 12)),
  
  # SI Fig 2
  tar_target(
    reed_muench_marginal_joint_bivariate_plot,
    plot_reed_muench_marginal_joint_bivariate(
      reed_muench_fit = ab_reed_muench_mcmc_antibody_mech
    )
  ),
  
  tar_target(file_reed_muench_marginal_joint_bivariate_plot,
             save_plot_pdf(reed_muench_marginal_joint_bivariate_plot,
                           "outputs/reed_muench/posterior_pairs_plot.pdf", 10, 10)),
  
  tar_target(
    phi_titer_summary_CVA6,
    make_phi_titer_summary(
      raw_df = ab_raw_CVA6,
      log_phi_draws = log_phi_draws_CVA6,
      virus = "CVA6"
    )
  ),
  
  tar_target(
    phi_titer_summary_EV71,
    make_phi_titer_summary(
      raw_df = ab_raw_EV71,
      log_phi_draws = log_phi_draws_EV71,
      virus = "EV71"
    )
  ),
  
  tar_target(
    phi_titer_summary_EV68,
    make_phi_titer_summary(
      raw_df = ab_raw_EV68,
      log_phi_draws = log_phi_draws_EV68,
      virus = "EV68"
    )
  ),
  
  tar_target(
    combined_phi_titer_summary,
    dplyr::bind_rows(
      phi_titer_summary_CVA6,
      phi_titer_summary_EV71,
      phi_titer_summary_EV68
    )
  ),
  
  tar_target(
    combined_phi_vs_titer_plot,
    plot_phi_vs_titer_combined(combined_phi_titer_summary)
  ),
  
  tar_target(file_combined_phi_vs_titer, save_plot_pdf(combined_phi_vs_titer_plot, "outputs/combined/combined_phi_vs_titer_plot.pdf", 10, 7), format = "file"),
  
  
  tar_target(
    conc_dilution_titer_data_CVA6,
    make_conc_dilution_titer_data(
      raw_df = ab_raw_CVA6,
      log_phi_draws = log_phi_draws_CVA6,
      virus = "CVA6"
    )
  ),
  
  tar_target(
    conc_dilution_titer_data_EV71,
    make_conc_dilution_titer_data(
      raw_df = ab_raw_EV71,
      log_phi_draws = log_phi_draws_EV71,
      virus = "EV71"
    )
  ),
  
  tar_target(
    conc_dilution_titer_data_EV68,
    make_conc_dilution_titer_data(
      raw_df = ab_raw_EV68,
      log_phi_draws = log_phi_draws_EV68,
      virus = "EV68"
    )
  ),
  
  tar_target(
    combined_conc_dilution_titer_data,
    dplyr::bind_rows(
      conc_dilution_titer_data_CVA6,
      conc_dilution_titer_data_EV71,
      conc_dilution_titer_data_EV68
    )
  ),
  
  tar_target(
    conc_vs_dilution_plot,
    plot_phi_over_dilution(combined_conc_dilution_titer_data)
  ),
  
  
  tar_target(
    fig_2,
    make_fig_2(
      reed_muench_plot,
      combined_phi_vs_titer_plot,
      conc_vs_dilution_plot
    )
  ),
  
  tar_target(file_fig_2, save_plot_pdf(fig_2, "outputs/combined/fig_2.pdf", 14, 6), format = "file"),
  
  tar_target(
    phi_age_group_data_CVA6,
    make_phi_age_group_data(
      serum_meta = ab_prep_CVA6$serum_meta,
      log_phi_draws = log_phi_draws_CVA6,
      virus = "CVA6"
    )
  ),
  
  tar_target(
    phi_age_group_data_EV71,
    make_phi_age_group_data(
      serum_meta = ab_prep_EV71$serum_meta,
      log_phi_draws = log_phi_draws_EV71,
      virus = "EV71"
    )
  ),
  
  tar_target(
    phi_age_group_data_EV68,
    make_phi_age_group_data(
      serum_meta = ab_prep_EV68$serum_meta,
      log_phi_draws = log_phi_draws_EV68,
      virus = "EV68"
    )
  ),
  
  tar_target(
    combined_phi_age_group_data,
    dplyr::bind_rows(
      phi_age_group_data_CVA6,
      phi_age_group_data_EV71,
      phi_age_group_data_EV68
    )
  ),
  
  tar_target(
    combined_phi_age_group_plot,
    plot_phi_by_age_group(combined_phi_age_group_data)
  ),
  
  tar_target(
    file_combined_phi_age_group,
    save_plot_pdf(
      combined_phi_age_group_plot,
      "outputs/combined/combined_phi_age_group_plot.pdf",
      12, 5
    ),
    format = "file"
  ),
  
  #### Seroreversion fits (testing) ####

  tar_target(
    sero_sr_data_CVA6,
    add_sero_switches(
      sero_stan_data_CVA6,
      estimate_omega = 1L,
      omega_prior_sd = sero_omega_prior_sd
    )
  ),
  
  tar_stan_mcmc( # SIS:
    name = "CVA6_sr",
    stan_files = "Stan/antibody_sero_sis.stan",
    data = sero_sr_data_CVA6,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  tar_stan_mcmc( # SIIS
    name = "CVA6_sr_3",
    stan_files = "Stan/antibody_sero_siis.stan",
    data = sero_sr_data_CVA6,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  tar_target(
    sero_sr_data_EV71,
    add_sero_switches(
      sero_stan_data_EV71,
      estimate_omega = 1L,
      omega_prior_sd = sero_omega_prior_sd
    )
  ),
  tar_stan_mcmc( # SIS:
    name = "EV71_sr",
    stan_files = "Stan/antibody_sero_sis.stan",
    data = sero_sr_data_EV71,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  tar_stan_mcmc( # SIIS
    name = "EV71_sr_3",
    stan_files = "Stan/antibody_sero_siis.stan",
    data = sero_sr_data_EV71,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  tar_target(
    sero_sr_data_EV68,
    add_sero_switches(
      sero_stan_data_EV68,
      estimate_omega = 1L,
      omega_prior_sd = sero_omega_prior_sd
    )
  ),
  
  tar_stan_mcmc( # SIS:
    name = "EV68_sr",
    stan_files = "Stan/antibody_sero_sis.stan",
    data = sero_sr_data_EV68,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  tar_stan_mcmc( # SIIS
    name = "EV68_sr_3",
    stan_files = "Stan/antibody_sero_siis.stan",
    data = sero_sr_data_EV68,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    seed = sero_seed,
    refresh = sero_refresh
  ),
  
  #### LOO-CV across SI/SII/SIS/SIIS models ####
  tar_target(
    CVA6_model_loo,
    compare_sero_models_loo(
      list(
        SI   = CVA6_two_comp_draws_antibody_sero_two_compartment,
        SII  = CVA6_draws_antibody_sero,
        SIS  = CVA6_sr_draws_antibody_sero_sis,
        SIIS = CVA6_sr_3_draws_antibody_sero_siis
      ),
      virus = "CVA6"
    )
  ),
  
  tar_target(
    EV71_model_loo,
    compare_sero_models_loo(
      list(
        SI   = EV71_two_comp_draws_antibody_sero_two_compartment,
        SII  = EV71_draws_antibody_sero,
        SIS  = EV71_sr_draws_antibody_sero_sis,
        SIIS = EV71_sr_3_draws_antibody_sero_siis
      ),
      virus = "EV71"
    )
  ),
  
  tar_target(
    EV68_model_loo,
    compare_sero_models_loo(
      list(
        SI   = EV68_draws_antibody_sero_two_compartment,   # note: no _two_comp prefix
        SII  = EV68_three_comp_draws_antibody_sero,
        SIS  = EV68_sr_draws_antibody_sero_sis,
        SIIS = EV68_sr_3_draws_antibody_sero_siis
      ),
      virus = "EV68"
    )
  ),
  
  tar_target(
    combined_model_loo_df,
    combine_model_loo_tables(CVA6_model_loo, EV71_model_loo, EV68_model_loo)
  ),
  
  #### k0 sensitivity analysis ####
  k0_targets,
  
  #### Seroreversion main results figures -> outputs/seroreversion ####
  
  tar_target(
    file_sr_sis_phi,
    save_sr_combined(
      "weighted_avg_phi", "SIS",
      sero_plot_data_CVA6, CVA6_sr_draws_antibody_sero_sis, CVA6_sr_summary_antibody_sero_sis,
      sero_plot_data_EV71, EV71_sr_draws_antibody_sero_sis, EV71_sr_summary_antibody_sero_sis,
      sero_plot_data_EV68, EV68_sr_draws_antibody_sero_sis, EV68_sr_summary_antibody_sero_sis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_siis_phi,
    save_sr_combined(
      "weighted_avg_phi", "SIIS",
      sero_plot_data_CVA6, CVA6_sr_3_draws_antibody_sero_siis, CVA6_sr_3_summary_antibody_sero_siis,
      sero_plot_data_EV71, EV71_sr_3_draws_antibody_sero_siis, EV71_sr_3_summary_antibody_sero_siis,
      sero_plot_data_EV68, EV68_sr_3_draws_antibody_sero_siis, EV68_sr_3_summary_antibody_sero_siis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_sis_sero,
    save_sr_combined(
      "serodynamics", "SIS",
      sero_plot_data_CVA6, CVA6_sr_draws_antibody_sero_sis, CVA6_sr_summary_antibody_sero_sis,
      sero_plot_data_EV71, EV71_sr_draws_antibody_sero_sis, EV71_sr_summary_antibody_sero_sis,
      sero_plot_data_EV68, EV68_sr_draws_antibody_sero_sis, EV68_sr_summary_antibody_sero_sis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_siis_sero,
    save_sr_combined(
      "serodynamics", "SIIS",
      sero_plot_data_CVA6, CVA6_sr_3_draws_antibody_sero_siis, CVA6_sr_3_summary_antibody_sero_siis,
      sero_plot_data_EV71, EV71_sr_3_draws_antibody_sero_siis, EV71_sr_3_summary_antibody_sero_siis,
      sero_plot_data_EV68, EV68_sr_3_draws_antibody_sero_siis, EV68_sr_3_summary_antibody_sero_siis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_sis_foi,
    save_sr_combined(
      "FOI_new_old", "SIS",
      sero_plot_data_CVA6, CVA6_sr_draws_antibody_sero_sis, CVA6_sr_summary_antibody_sero_sis,
      sero_plot_data_EV71, EV71_sr_draws_antibody_sero_sis, EV71_sr_summary_antibody_sero_sis,
      sero_plot_data_EV68, EV68_sr_draws_antibody_sero_sis, EV68_sr_summary_antibody_sero_sis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_siis_foi,
    save_sr_combined(
      "FOI_new_old", "SIIS",
      sero_plot_data_CVA6, CVA6_sr_3_draws_antibody_sero_siis, CVA6_sr_3_summary_antibody_sero_siis,
      sero_plot_data_EV71, EV71_sr_3_draws_antibody_sero_siis, EV71_sr_3_summary_antibody_sero_siis,
      sero_plot_data_EV68, EV68_sr_3_draws_antibody_sero_siis, EV68_sr_3_summary_antibody_sero_siis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_sis_omega,
    save_sr_combined(
      "omega", "SIS",
      sero_plot_data_CVA6, CVA6_sr_draws_antibody_sero_sis, CVA6_sr_summary_antibody_sero_sis,
      sero_plot_data_EV71, EV71_sr_draws_antibody_sero_sis, EV71_sr_summary_antibody_sero_sis,
      sero_plot_data_EV68, EV68_sr_draws_antibody_sero_sis, EV68_sr_summary_antibody_sero_sis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_siis_omega,
    save_sr_combined(
      "omega", "SIIS",
      sero_plot_data_CVA6, CVA6_sr_3_draws_antibody_sero_siis, CVA6_sr_3_summary_antibody_sero_siis,
      sero_plot_data_EV71, EV71_sr_3_draws_antibody_sero_siis, EV71_sr_3_summary_antibody_sero_siis,
      sero_plot_data_EV68, EV68_sr_3_draws_antibody_sero_siis, EV68_sr_3_summary_antibody_sero_siis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_param_tables,
    save_sr_param_tables(list(
      list(virus = "CVA6", model = "SIS",  summary = CVA6_sr_summary_antibody_sero_sis),
      list(virus = "CVA6", model = "SIIS", summary = CVA6_sr_3_summary_antibody_sero_siis),
      list(virus = "EV71", model = "SIS",  summary = EV71_sr_summary_antibody_sero_sis),
      list(virus = "EV71", model = "SIIS", summary = EV71_sr_3_summary_antibody_sero_siis),
      list(virus = "EV68", model = "SIS",  summary = EV68_sr_summary_antibody_sero_sis),
      list(virus = "EV68", model = "SIIS", summary = EV68_sr_3_summary_antibody_sero_siis)
    )),
    format = "file"
  ),
  
  #### vectorised sigma ####
  tar_target(sero_stan_data_CVA6_shared,
             make_sero_stan_data(serum_meta_CVA6, log_phi_draws_CVA6,
                                 sigma_by_state = 0L)),
  
  tar_stan_mcmc(
    name = "CVA6_shared",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_stan_data_CVA6_shared,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh
  ),
  
  tar_target(
    sigma_loo_CVA6,
    compare_sero_models_loo(list(
      shared_sigma   = CVA6_shared_draws_antibody_sero,
      by_state_sigma = CVA6_draws_antibody_sero
    ), virus = "CVA6")
  ),
  
  tar_target(
    sigma_reclass_CVA6,
    serostate_reclassification(CVA6_shared_summary_antibody_sero,
                               CVA6_summary_antibody_sero,
                               n = sero_plot_data_CVA6$n, k = 3)
  ),
  
  # EV68
  tar_target(sero_stan_data_EV68_shared,
             make_sero_stan_data(serum_meta_EV68, log_phi_draws_EV68,
                                 sigma_by_state = 0L)),
  
  tar_stan_mcmc(
    name = "EV68_shared",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_stan_data_EV68_shared,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh
  ),
  
  tar_target(
    sigma_loo_EV68,
    compare_sero_models_loo(list(
      shared_sigma   = EV68_shared_draws_antibody_sero,
      by_state_sigma = EV68_draws_antibody_sero
    ), virus = "EV68")
  ),
  
  tar_target(
    sigma_reclass_EV68,
    serostate_reclassification(EV68_shared_summary_antibody_sero,
                               EV68_summary_antibody_sero,
                               n = sero_plot_data_EV68$n, k = 3)
  ),
  
  # EV71
  tar_target(sero_stan_data_EV71_shared,
             make_sero_stan_data(serum_meta_EV71, log_phi_draws_EV71,
                                 sigma_by_state = 0L)),
  
  tar_stan_mcmc(
    name = "EV71_shared",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_stan_data_EV71_shared,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh
  ),
  
  tar_target(
    sigma_loo_EV71,
    compare_sero_models_loo(list(
      shared_sigma   = EV71_shared_draws_antibody_sero,
      by_state_sigma = EV71_draws_antibody_sero
    ), virus = "EV71")
  ),
  
  tar_target(
    sigma_reclass_EV71,
    serostate_reclassification(EV71_shared_summary_antibody_sero,
                               EV71_summary_antibody_sero,
                               n = sero_plot_data_EV71$n, k = 3)
  )
)
