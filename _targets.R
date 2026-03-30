# _targets.R
library(targets)
library(stantargets)

# global options
tar_option_set(
  packages = c(
    "targets", "stantargets", "cmdstanr",
    "dplyr", "tidyr", "purrr", "stringr",
    "posterior", "ggplot2", "gt",
    "matrixStats", "reshape2", "cowplot", "gridExtra", "grid"
  ),
  format = "qs"
)

# load all functions in R/ directory
tar_source()

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

# thin antibody posterior draws for speed
n_phi_draws_thin <- 500


## (need to figure out how to set this via relative paths)
# # set device path to Google Chrome for gt tables 
# Sys.setenv(CHROMOTE_CHROME = "C:/Users/alrefae/AppData/Local/Google/Chrome/Application/chrome.exe")


# ---------- Pipeline ----------
list(
  # -------- Antibody data prep (raw well observations -> stan_data + serum_meta) --------
  tar_target(ab_prep_EV71, prep_antibody_data(virus = "EV71", raw_prefix = "EV71")),
  tar_target(ab_prep_CVA6, prep_antibody_data(virus = "CVA6", raw_prefix = "CA6")),
  tar_target(ab_prep_EV68, prep_antibody_data(virus = "EV68", raw_prefix = "E68")),
  
  # -------- Antibody fits (mechanistic model) --------
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
  tar_target(ab_raw_EV71, read.csv("data/EV71_raw_well_observations.csv")),
  tar_target(ab_raw_CVA6, read.csv("data/CA6_raw_well_observations.csv")),
  tar_target(ab_raw_EV68, read.csv("data/E68_raw_well_observations.csv")),
  
  # --- Prior + Posterior PPC overlay plots ---
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
  
  # --- Save PPC plots to file ---
  tar_target(
    file_ab_ppc_EV71,
    {
      save_plot_pdf(ab_ppc_EV71, "outputs/EV71/EV71_antibody_ppc_overlay.pdf", width = 12, height = 8)
    },
    format = "file"
  ),
  
  tar_target(
    file_ab_ppc_CVA6,
    {
      save_plot_pdf(ab_ppc_CVA6, "outputs/CVA6/CVA6_antibody_ppc_overlay.pdf", width = 12, height = 8)
    },
    format = "file"
  ),
  
  tar_target(
    file_ab_ppc_EV68,
    {
      save_plot_pdf(ab_ppc_EV68, "outputs/EV68/EV68_antibody_ppc_overlay.pdf", width = 12, height = 8)
    },
    format = "file"
  ),
  
  # -------- Extract full posterior draws of log(phi) per individual (thinning) --------
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
    log_phi_draws_CVA6,
    extract_log_phi_draws(
      draws_obj = ab_CVA6_draws_antibody_mech,
      n_individuals = ab_prep_CVA6$stan_data$n_individuals,
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
    sero_stan_data_EV71,
    make_sero_stan_data(
      serum_meta = ab_prep_EV71$serum_meta,
      log_phi_draws = log_phi_draws_EV71
    )
  ),
  
  tar_target(
    sero_stan_data_CVA6,
    make_sero_stan_data(
      serum_meta = ab_prep_CVA6$serum_meta,
      log_phi_draws = log_phi_draws_CVA6
    )
  ),
  
  tar_target(
    sero_stan_data_EV68,
    make_sero_stan_data(
      serum_meta = ab_prep_EV68$serum_meta,
      log_phi_draws = log_phi_draws_EV68
    )
  ),
  
  # -------- Build plotting data --------
  tar_target(
    sero_plot_data_EV71,
    make_sero_plot_data(
      serum_meta = ab_prep_EV71$serum_meta,
      log_phi_draws = log_phi_draws_EV71
    )
  ),
  
  tar_target(
    sero_plot_data_CVA6,
    make_sero_plot_data(
      serum_meta = ab_prep_CVA6$serum_meta,
      log_phi_draws = log_phi_draws_CVA6
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
  tar_stan_mcmc(
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
  
  tar_stan_mcmc(
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
  
  tar_stan_mcmc(
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
  
  tar_stan_mcmc(
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
      EV71_draws_antibody_sero,
      CVA6_draws_antibody_sero,
      EV68_draws_antibody_sero_two_compartment
    )
  ),
  
  tar_target(
    summaries_combined,
    combine_summaries_three_viruses(
      EV71_summary_antibody_sero,
      CVA6_summary_antibody_sero,
      EV68_summary_antibody_sero_two_compartment
    )
  ),
  
  # -------- FOI plots --------
  tar_target(
    EV71_FOI_plot,
    plot_foi_short_long("EV71", sero_plot_data_EV71, EV71_summary_antibody_sero,
                        EV71_draws_antibody_sero)
  ),
  tar_target(
    CVA6_FOI_plot,
    plot_foi_short_long("CVA6", sero_plot_data_CVA6, CVA6_summary_antibody_sero,
                        CVA6_draws_antibody_sero)
  ),
  tar_target(
    EV68_FOI_plot,
    plot_foi_short_long("EV68", sero_plot_data_EV68, EV68_summary_antibody_sero_two_compartment,
                        EV68_draws_antibody_sero_two_compartment)
  ),
  
  tar_target(
    combined_FOI_plot,
    plot_foi_combined(draws_combined, max_age = max(sero_plot_data_EV71$age_max, sero_plot_data_CVA6$age_max))
  ),
  
  # -------- Posterior density plots (combined) --------
  tar_target(
    combined_density_plots,
    plot_param_densities_combined(draws_combined)
  ),
  
  # -------- Parameter tables --------
  tar_target(
    EV71_param_table,
    make_param_table("EV71", EV71_summary_antibody_sero)
  ),
  tar_target(
    CVA6_param_table,
    make_param_table("CVA6", CVA6_summary_antibody_sero)
  ),
  tar_target(
    EV68_param_table,
    make_param_table("EV68", EV68_summary_antibody_sero_two_compartment)
  ),
  
  # -------- Phi plots --------
  tar_target(
    EV71_phi_plot,
    plot_phi_points("EV71", sero_plot_data_EV71, EV71_summary_antibody_sero)
  ),
  tar_target(
    CVA6_phi_plot,
    plot_phi_points("CVA6", sero_plot_data_CVA6, CVA6_summary_antibody_sero)
  ),
  tar_target(
    EV68_phi_plot,
    plot_phi_points("EV68", sero_plot_data_EV68, EV68_summary_antibody_sero_two_compartment)
  ),
  
  # -------- Serodynamics plots --------
  tar_target(
    EV71_serodynamics_plot,
    plot_serodynamics("EV71", sero_plot_data_EV71, EV71_draws_antibody_sero)
  ),
  tar_target(
    CVA6_serodynamics_plot,
    plot_serodynamics("CVA6", sero_plot_data_CVA6, CVA6_draws_antibody_sero)
  ),
  tar_target(
    EV68_serodynamics_plot,
    plot_serodynamics("EV68", sero_plot_data_EV68, EV68_draws_antibody_sero_two_compartment)
  ),
  
  # -------- Weighted-average phi plot --------
  tar_target(
    EV71_weighted_avg_phi_plot,
    plot_weighted_avg_phi("EV71", sero_plot_data_EV71,
                          EV71_draws_antibody_sero,
                          EV71_summary_antibody_sero)
  ),
  tar_target(
    CVA6_weighted_avg_phi_plot,
    plot_weighted_avg_phi("CVA6", sero_plot_data_CVA6,
                          CVA6_draws_antibody_sero,
                          CVA6_summary_antibody_sero)
  ),
  tar_target(
    EV68_weighted_avg_phi_plot,
    plot_weighted_avg_phi("EV68", sero_plot_data_EV68,
                          EV68_draws_antibody_sero_two_compartment,
                          EV68_summary_antibody_sero_two_compartment)
  ),
  
  # -------- Save outputs --------
  tar_target(file_EV71_FOI, save_plot_pdf(EV71_FOI_plot, "outputs/EV71/EV71_FOI_plot.pdf", 8, 6), format = "file"),
  tar_target(file_CVA6_FOI, save_plot_pdf(CVA6_FOI_plot, "outputs/CVA6/CVA6_FOI_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_FOI, save_plot_pdf(EV68_FOI_plot, "outputs/EV68/EV68_FOI_plot.pdf", 8, 6), format = "file"),
  tar_target(file_combined_FOI, save_plot_pdf(combined_FOI_plot, "outputs/combined/combined_FOI_plot.pdf", 8, 6), format = "file"),
  
  tar_target(file_density, save_plot_pdf(combined_density_plots, "outputs/combined/combined_density_plots.pdf", 10, 7), format = "file"),
  
  tar_target(file_EV71_phi, save_plot_pdf(EV71_phi_plot, "outputs/EV71/EV71_phi_plot.pdf", 8, 6), format = "file"),
  tar_target(file_CVA6_phi, save_plot_pdf(CVA6_phi_plot, "outputs/CVA6/CVA6_phi_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_phi, save_plot_pdf(EV68_phi_plot, "outputs/EV68/EV68_phi_plot.pdf", 8, 6), format = "file"),
  
  tar_target(file_EV71_serodynamics, save_plot_pdf(EV71_serodynamics_plot, "outputs/EV71/EV71_serodynamics_plot.pdf", 8, 6), format = "file"),
  tar_target(file_CVA6_serodynamics, save_plot_pdf(CVA6_serodynamics_plot, "outputs/CVA6/CVA6_serodynamics_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_serodynamics, save_plot_pdf(EV68_serodynamics_plot, "outputs/EV68/EV68_serodynamics_plot.pdf", 8, 6), format = "file"),
  
  tar_target(file_EV71_weighted, save_plot_pdf(EV71_weighted_avg_phi_plot, "outputs/EV71/EV71_weighted_avg_phi_plot.pdf", 10, 7), format = "file"),
  tar_target(file_CVA6_weighted, save_plot_pdf(CVA6_weighted_avg_phi_plot, "outputs/CVA6/CVA6_weighted_avg_phi_plot.pdf", 10, 7), format = "file"),
  tar_target(file_EV68_weighted, save_plot_pdf(EV68_weighted_avg_phi_plot, "outputs/EV68/EV68_weighted_avg_phi_plot.pdf", 10, 7), format = "file"),
  
  tar_target(file_EV71_table, save_gt_pdf(EV71_param_table, "outputs/EV71/EV71_param_table.pdf"), format = "file"),
  tar_target(file_CVA6_table, save_gt_pdf(CVA6_param_table, "outputs/CVA6/CVA6_param_table.pdf"), format = "file"),
  tar_target(file_EV68_table, save_gt_pdf(EV68_param_table, "outputs/EV68/EV68_param_table.pdf"), format = "file")
)
