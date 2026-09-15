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
    "priorsense", "tibble",
    
    # plotting / saving pdfs
    "ggplot2", "gt", "gtExtras", "cowplot", "grid", "gridExtra", "khroma",
    "patchwork", "scales", "pdftools", "showtext", "viridis", "GGally"
  ),
  format = "qs",
  resources = tar_resources(
      qs = tar_resources_qs(preset = "high")
    ),
  memory = "transient",
  garbage_collection = TRUE
)

# ---------- Global MCMC settings ----------

# antibody model
ab_iter_warmup   <- 500 #5
ab_iter_sampling <- 500 #5
ab_chains        <- 4 #1
ab_parallel      <- 4 #1
ab_seed          <- 1
ab_refresh       <- 200 #1

# serocatalytic model
sero_iter_sampling <- 1000 #2
sero_iter_warmup   <- floor(sero_iter_sampling / 5) #1
sero_chains        <- 4 #1
sero_parallel      <- 4 #1
sero_seed          <- 2
sero_refresh       <- 200 #1

# seroreversion fits: omega is weakly identified and needs longer adaptation
sero_sr_iter_warmup <- 800 #2
sero_sr_adapt_delta <- 0.9
sero_omega_prior_sd <- 0.5
sero_sigma_by_state <- 1L

# thin antibody posterior draws for speed (NULL = no thinning)
n_phi_draws_thin <- 250 #NULL 

# k0/k1/phi sensitivity analysis
K0_SPECS <- c("-2", "-1", "0", "1", "2", "free")
K0_IDS   <- c("m2", "m1", "z0", "p1", "p2", "free")
K0_FIXED <- c("-2", "-1", "0", "1", "2")
K0_PHI_PRIOR_SCALE <- 500
K0_PRIOR_SD        <- 5

K0_S2_THIN     <- n_phi_draws_thin
K0_S2_CHAINS   <- sero_chains
K0_S2_WARMUP   <- 200
K0_S2_SAMPLING <- 200
K0_S2_REFRESH  <- 20

K1_PRIOR_SCALE <- 1

# threading
n_cores      <- parallel::detectCores(logical = FALSE)
sero_threads <- max(1L, floor(n_cores / sero_chains))
ab_threads   <- max(1L, floor(n_cores / ab_chains))
stan_cpp     <- list(stan_threads = TRUE)

sero_vars <- list(
  sii = c("mu", "mu_x", "sigma_phi", "sigma_x", "psi",
          "lambda_1", "kappa", "lambda_long",
          "log_likelihood", "prob_by_group", "lp__"),
  
  si  = c("mu", "mu_x", "sigma_phi", "sigma_x",
          "lambda_1", "kappa", "lambda_long",
          "log_likelihood", "prob_by_group", "lp__"),
  
  siis = c("mu", "mu_x", "sigma_phi", "sig_x", "psi", "omega",
           "lambda_1", "kappa", "lambda_long", "state_probs",
           "mean_years_to_reversion", "lprior",
           "log_lik", "prob_by_group", "lp__"),
  
  sis  = c("mu", "mu_x", "sigma_phi", "sig_x", "omega",
           "lambda_1", "kappa", "lambda_long", "state_probs",
           "mean_years_to_reversion", "lprior",
           "log_lik", "prob_by_group", "lp__")
)

# load all functions in R/ directory
tar_source()
source("k0_pipeline.R")
source("phi_k1_prior_pipeline.R")

# compile .stan files
invisible(lapply(
  c("Stan/antibody_mech.stan", "Stan/antibody_mech_k0.stan",
    #"Stan/antibody_mech_phi_k1.stan",
    "Stan/antibody_sero.stan", "Stan/antibody_sero_two_compartment.stan",
    "Stan/antibody_sero_sis.stan", "Stan/antibody_sero_siis.stan"),
  function(f) cmdstanr::cmdstan_model(f, cpp_options = stan_cpp)
))


# to render tables
ensure_chromote_browser()

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
    refresh = ab_refresh,
    
    cpp_options = stan_cpp,
    threads_per_chain = ab_threads
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
    refresh = ab_refresh,
    
    cpp_options = stan_cpp,
    threads_per_chain = ab_threads
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
    refresh = ab_refresh,
    
    cpp_options = stan_cpp,
    threads_per_chain = ab_threads
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
    refresh = sero_refresh,
    
    variables = sero_vars$sii,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
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
    refresh = sero_refresh,
    
    variables = sero_vars$sii,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
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
    refresh = sero_refresh,
    
    variables = sero_vars$si,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
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
    refresh = sero_refresh,
    
    variables = sero_vars$sii,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
  ),
  
  # -------- Combine draws & summaries --------
  tar_target(
    draws_combined,
    combine_draws_three_viruses(
      CVA6_shared_draws_antibody_sero,
      EV71_shared_draws_antibody_sero,
      EV68_draws_antibody_sero_two_compartment
    )
  ),
  
  tar_target(
    summaries_combined,
    combine_summaries_three_viruses(
      CVA6_shared_summary_antibody_sero,
      EV71_shared_summary_antibody_sero,
      EV68_summary_antibody_sero_two_compartment
    )
  ),
  
  # -------- Phi plots --------
  tar_target(
    CVA6_phi_plot,
    plot_phi_points("CVA6", sero_plot_data_CVA6, CVA6_shared_summary_antibody_sero)
  ),
  tar_target(
    EV71_phi_plot,
    plot_phi_points("EV71", sero_plot_data_EV71, EV71_shared_summary_antibody_sero)
  ),
  tar_target(
    EV68_phi_plot,
    plot_phi_points("EV68", sero_plot_data_EV68, EV68_summary_antibody_sero_two_compartment)
  ),
  
  # -------- Weighted-average phi plot --------
  tar_target(
    CVA6_weighted_avg_phi_plot,
    plot_weighted_avg_phi("CVA6", sero_plot_data_CVA6,
                          CVA6_shared_draws_antibody_sero,
                          CVA6_shared_summary_antibody_sero)
  ),
  tar_target(
    EV71_weighted_avg_phi_plot,
    plot_weighted_avg_phi("EV71", sero_plot_data_EV71,
                          EV71_shared_draws_antibody_sero,
                          EV71_shared_summary_antibody_sero)
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
    plot_serodynamics("CVA6", sero_plot_data_CVA6, CVA6_shared_draws_antibody_sero)
  ),
  tar_target(
    EV71_serodynamics_plot,
    plot_serodynamics("EV71", sero_plot_data_EV71, EV71_shared_draws_antibody_sero)
  ),
  tar_target( # SI model (serostate-specific sigma_phi)
    EV68_serodynamics_plot,
    plot_serodynamics("EV68", sero_plot_data_EV68, EV68_draws_antibody_sero_two_compartment, three_comp=FALSE)
  ),
  
  tar_target( # SII model (shared sigma_phi)
    EV68_serodynamics_plot_three_comp,
    plot_serodynamics("EV68", sero_plot_data_EV68, EV68_three_comp_shared_draws_antibody_sero)
  ),
  
  tar_target( # EV-D68 SI stratified by year
    EV68_serodynamics_by_year_plot,
    plot_serodynamics_by_year(
      "EV68", sero_plot_data_EV68,
      strata_draws = list(
        `2006` = EV68_yr2006_draws_antibody_sero_two_compartment,
        `2011` = EV68_yr2011_draws_antibody_sero_two_compartment,
        `2017` = EV68_yr2017_draws_antibody_sero_two_compartment))
  ),
  
  # -------- FOI plots --------
  tar_target(
    CVA6_FOI_plot,
    plot_foi_short_long("CVA6", sero_plot_data_CVA6, CVA6_shared_draws_antibody_sero)
  ),
  tar_target(
    EV71_FOI_plot,
    plot_foi_short_long("EV71", sero_plot_data_EV71, EV71_shared_draws_antibody_sero)
  ),
  tar_target(
    EV68_FOI_plot,
    plot_foi_short_long("EV68", sero_plot_data_EV68, EV68_draws_antibody_sero_two_compartment)
  ),
  
  tar_target( # stratified by survey year
    EV68_foi_by_year_plot,
    plot_foi_by_year("EV68", sero_plot_data_EV68,
                     strata_draws = list(
                       `2006` = EV68_yr2006_draws_antibody_sero_two_compartment,
                       `2011` = EV68_yr2011_draws_antibody_sero_two_compartment,
                       `2017` = EV68_yr2017_draws_antibody_sero_two_compartment))
  ),
  
  # single plot
  tar_target(
    combined_FOI_plot,
    plot_foi_combined(draws_combined, max_age = max(sero_plot_data_CVA6$age_max,
                                                    sero_plot_data_EV71$age_max,
                                                    sero_plot_data_EV68$age_max))
  ),
  
  
  # new vs old FOI plots
  
  tar_target(old_foi_CVA6, summarise_old_foi("CVA6", sero_plot_data_CVA6$age_max)),
  tar_target(old_foi_EV71, summarise_old_foi("EV71", sero_plot_data_EV71$age_max)),
  tar_target(old_foi_EV68, summarise_old_foi("EV68", sero_plot_data_EV68$age_max)),
  
  # -------- Parameter tables --------
  tar_target(
    CVA6_param_table,
    make_param_table("CVA6", CVA6_shared_summary_antibody_sero, CVA6_shared_draws_antibody_sero)
  ),
  tar_target(
    EV71_param_table,
    make_param_table("EV71", EV71_shared_summary_antibody_sero, EV71_shared_draws_antibody_sero)
  ),
  tar_target(
    EV68_param_table,
    make_param_table("EV68", EV68_summary_antibody_sero_two_compartment, EV68_draws_antibody_sero_two_compartment)
  ),
  
  # -------- Save outputs --------
  tar_target(file_ab_ppc_CVA6, save_plot_pdf(ab_ppc_CVA6, "outputs/CVA6/CVA6_antibody_ppc_overlay.pdf", width = 12, height = 8), format = "file"),
  tar_target(file_ab_ppc_EV71, save_plot_pdf(ab_ppc_EV71, "outputs/EV71/EV71_antibody_ppc_overlay.pdf", width = 12, height = 8), format = "file"),
  tar_target(file_ab_ppc_EV68, save_plot_pdf(ab_ppc_EV68, "outputs/EV68/EV68_antibody_ppc_overlay.pdf", width = 12, height = 8), format = "file"),
  
  tar_target(file_CVA6_FOI, save_plot_pdf(CVA6_FOI_plot, "outputs/CVA6/CVA6_FOI_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV71_FOI, save_plot_pdf(EV71_FOI_plot, "outputs/EV71/EV71_FOI_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_FOI, save_plot_pdf(EV68_FOI_plot, "outputs/EV68/EV68_FOI_plot.pdf", 8, 6), format = "file"),
  
  tar_target(file_combined_foi, {
    g <- plot_patchwork(
      plot_foi_new_old("CVA6", foi_new_summary(CVA6_shared_draws_antibody_sero,
                                               sero_plot_data_CVA6$age_max), old_foi_CVA6),
      plot_foi_new_old("EV71", foi_new_summary(EV71_shared_draws_antibody_sero,
                                               sero_plot_data_EV71$age_max), old_foi_EV71),
      plot_foi_new_old("EV68", foi_new_summary(EV68_draws_antibody_sero_two_compartment,
                                               sero_plot_data_EV68$age_max), old_foi_EV68),
      FOI = TRUE,
      fourth_plot = EV68_foi_by_year_plot
    )
      save_plot_pdf(g, "outputs/combined/combined_FOI_grid_plot.pdf", 12, 10)
    }, format = "file"
  ),
  
  tar_target(
    file_density,
    save_plot_pdf(
      plot_param_densities_combined(draws_combined),
      "outputs/combined/combined_density_plots.pdf", 12, 10
    ),
    format = "file"
  ),
  
  tar_target(file_CVA6_phi, save_plot_pdf(CVA6_phi_plot, "outputs/CVA6/CVA6_phi_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV71_phi, save_plot_pdf(EV71_phi_plot, "outputs/EV71/EV71_phi_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_phi, save_plot_pdf(EV68_phi_plot, "outputs/EV68/EV68_phi_plot.pdf", 8, 6), format = "file"),
  
  tar_target(file_CVA6_serodynamics, save_plot_pdf(CVA6_serodynamics_plot, "outputs/CVA6/CVA6_serodynamics_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV71_serodynamics, save_plot_pdf(EV71_serodynamics_plot, "outputs/EV71/EV71_serodynamics_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_serodynamics, save_plot_pdf(EV68_serodynamics_plot, "outputs/EV68/EV68_serodynamics_plot.pdf", 8, 6), format = "file"),
  tar_target(file_EV68_serodynamics_three_comp, save_plot_pdf(EV68_serodynamics_plot_three_comp, "outputs/EV68/EV68_serodynamics_plot_three_comp.pdf", 8, 6), format = "file"),
  
  tar_target(
    file_combined_serodynamics,
    {
      g <- plot_patchwork(
        CVA6_serodynamics_plot,
        EV71_serodynamics_plot,
        EV68_serodynamics_plot,
        lines=FALSE,
        fourth_plot = EV68_serodynamics_by_year_plot
      )
      save_plot_pdf(g, "outputs/combined/combined_serodynamics_plot.pdf", 12, 10)
    },
    format = "file"
  ),
  
  tar_target(file_CVA6_weighted, save_plot_pdf(CVA6_weighted_avg_phi_plot, "outputs/CVA6/CVA6_weighted_avg_phi_plot.pdf", 10, 7), format = "file"),
  tar_target(file_EV71_weighted, save_plot_pdf(EV71_weighted_avg_phi_plot, "outputs/EV71/EV71_weighted_avg_phi_plot.pdf", 10, 7), format = "file"),
  tar_target(file_EV68_weighted, save_plot_pdf(EV68_weighted_avg_phi_plot, "outputs/EV68/EV68_weighted_avg_phi_plot.pdf", 10, 7), format = "file"),
  
  tar_target(
    file_combined_weighted_avg_phi,
    {
      g <- plot_patchwork(
        CVA6_weighted_avg_phi_plot,
        EV71_weighted_avg_phi_plot,
        EV68_weighted_avg_phi_plot
      )
      save_plot_pdf(g, "outputs/combined/combined_weighted_avg_phi_plot.pdf", 12, 10)
    },
    format = "file"
  ),
  
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
    refresh = 200,
    
    cpp_options = stan_cpp,
    threads_per_chain = 1L
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

  tar_stan_mcmc(
    name = "CVA6_two_comp",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_stan_data_CVA6,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains,
    parallel_chains = sero_parallel,
    seed = sero_seed,
    refresh = sero_refresh,

    variables = sero_vars$si,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
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
    refresh = sero_refresh,

    variables = sero_vars$si,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
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
  
  tar_target(file_survival_plots, save_plot_pdf(survival_plots, "outputs/misc/survival_plots.pdf", 8, 12), format="file"),
  
  # SI Fig 2
  tar_target(
    reed_muench_marginal_joint_bivariate_plot,
    plot_reed_muench_marginal_joint_bivariate(
      reed_muench_fit = ab_reed_muench_mcmc_antibody_mech
    )
  ),
  
  tar_target(file_reed_muench_marginal_joint_bivariate_plot,
             save_plot_pdf(reed_muench_marginal_joint_bivariate_plot,
                           "outputs/reed_muench/posterior_pairs_plot.pdf", 10, 10), format="file"),
  
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
  
  #### Seroreversion fits ####

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
    refresh = sero_refresh,
    
    variables = sero_vars$sis,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
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
    refresh = sero_refresh,
    
    variables = sero_vars$siis,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
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
    refresh = sero_refresh,
    
    variables = sero_vars$sis,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
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
    refresh = sero_refresh,
    
    variables = sero_vars$siis,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
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
    refresh = sero_refresh,
    
    variables = sero_vars$sis,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
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
    refresh = sero_refresh,
    
    variables = sero_vars$siis,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
  ),
  
  #### k0 sensitivity analysis ####
  
  k0_targets,
  
  #### k1 / phi prior sensitivity analysis ####
  ps_targets,
  
  #### Seroreversion main results figures -> outputs/seroreversion ####
  
  tar_target(
    file_sr_sis_phi,
    save_sr_combined(
      "weighted_avg_phi", "SIS",
      sero_plot_data_CVA6, CVA6_sr_shared_draws_antibody_sero_sis, CVA6_sr_shared_summary_antibody_sero_sis,
      sero_plot_data_EV71, EV71_sr_shared_draws_antibody_sero_sis, EV71_sr_shared_summary_antibody_sero_sis,
      sero_plot_data_EV68, EV68_sr_shared_draws_antibody_sero_sis, EV68_sr_shared_summary_antibody_sero_sis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_siis_phi,
    save_sr_combined(
      "weighted_avg_phi", "SIIS",
      sero_plot_data_CVA6, CVA6_sr_3_shared_draws_antibody_sero_siis, CVA6_sr_3_shared_summary_antibody_sero_siis,
      sero_plot_data_EV71, EV71_sr_3_shared_draws_antibody_sero_siis, EV71_sr_3_shared_summary_antibody_sero_siis,
      sero_plot_data_EV68, EV68_sr_3_shared_draws_antibody_sero_siis, EV68_sr_3_shared_summary_antibody_sero_siis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_sis_sero,
    save_sr_combined(
      "serodynamics", "SIS",
      sero_plot_data_CVA6, CVA6_sr_shared_draws_antibody_sero_sis, CVA6_sr_shared_summary_antibody_sero_sis,
      sero_plot_data_EV71, EV71_sr_shared_draws_antibody_sero_sis, EV71_sr_shared_summary_antibody_sero_sis,
      sero_plot_data_EV68, EV68_sr_shared_draws_antibody_sero_sis, EV68_sr_shared_summary_antibody_sero_sis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_siis_sero,
    save_sr_combined(
      "serodynamics", "SIIS",
      sero_plot_data_CVA6, CVA6_sr_3_shared_draws_antibody_sero_siis, CVA6_sr_3_shared_summary_antibody_sero_siis,
      sero_plot_data_EV71, EV71_sr_3_shared_draws_antibody_sero_siis, EV71_sr_3_shared_summary_antibody_sero_siis,
      sero_plot_data_EV68, EV68_sr_3_shared_draws_antibody_sero_siis, EV68_sr_3_shared_summary_antibody_sero_siis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_sis_foi,
    save_sr_combined(
      "FOI_new_old", "SIS",
      sero_plot_data_CVA6, CVA6_sr_shared_draws_antibody_sero_sis, CVA6_sr_shared_summary_antibody_sero_sis,
      sero_plot_data_EV71, EV71_sr_shared_draws_antibody_sero_sis, EV71_sr_shared_summary_antibody_sero_sis,
      sero_plot_data_EV68, EV68_sr_shared_draws_antibody_sero_sis, EV68_sr_shared_summary_antibody_sero_sis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_siis_foi,
    save_sr_combined(
      "FOI_new_old", "SIIS",
      sero_plot_data_CVA6, CVA6_sr_3_shared_draws_antibody_sero_siis, CVA6_sr_3_shared_summary_antibody_sero_siis,
      sero_plot_data_EV71, EV71_sr_3_shared_draws_antibody_sero_siis, EV71_sr_3_shared_summary_antibody_sero_siis,
      sero_plot_data_EV68, EV68_sr_3_shared_draws_antibody_sero_siis, EV68_sr_3_shared_summary_antibody_sero_siis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_sis_omega,
    save_sr_combined(
      "omega", "SIS",
      sero_plot_data_CVA6, CVA6_sr_shared_draws_antibody_sero_sis, CVA6_sr_shared_summary_antibody_sero_sis,
      sero_plot_data_EV71, EV71_sr_shared_draws_antibody_sero_sis, EV71_sr_shared_summary_antibody_sero_sis,
      sero_plot_data_EV68, EV68_sr_shared_draws_antibody_sero_sis, EV68_sr_shared_summary_antibody_sero_sis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_siis_omega,
    save_sr_combined(
      "omega", "SIIS",
      sero_plot_data_CVA6, CVA6_sr_3_shared_draws_antibody_sero_siis, CVA6_sr_3_shared_summary_antibody_sero_siis,
      sero_plot_data_EV71, EV71_sr_3_shared_draws_antibody_sero_siis, EV71_sr_3_shared_summary_antibody_sero_siis,
      sero_plot_data_EV68, EV68_sr_3_shared_draws_antibody_sero_siis, EV68_sr_3_shared_summary_antibody_sero_siis,
      omega_prior_sd = sero_omega_prior_sd
    ),
    format = "file"
  ),
  
  tar_target(
    file_sr_param_tables,
    save_sr_param_tables(list(
      list(virus = "CVA6", model = "SIS",  summary = CVA6_sr_shared_summary_antibody_sero_sis),
      list(virus = "CVA6", model = "SIIS", summary = CVA6_sr_3_shared_summary_antibody_sero_siis),
      list(virus = "EV71", model = "SIS",  summary = EV71_sr_shared_summary_antibody_sero_sis),
      list(virus = "EV71", model = "SIIS", summary = EV71_sr_3_shared_summary_antibody_sero_siis),
      list(virus = "EV68", model = "SIS",  summary = EV68_sr_shared_summary_antibody_sero_sis),
      list(virus = "EV68", model = "SIIS", summary = EV68_sr_3_shared_summary_antibody_sero_siis)
    )),
    format = "file"
  ),
  
  #### vectorised sigma ####
  
  # CVA6
  tar_target(sero_stan_data_CVA6_shared,
             make_sero_stan_data(ab_prep_CVA6$serum_meta, log_phi_draws_CVA6,
                                 sigma_by_state = 0L)),
  
  # SII
  tar_stan_mcmc(
    name = "CVA6_shared",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_stan_data_CVA6_shared,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$sii,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
  ),
  
  # SI
  tar_stan_mcmc(
    name = "CVA6_two_comp_shared",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_stan_data_CVA6_shared,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$si,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_target(sero_sr_shared_data_CVA6,
             add_sero_switches(sero_stan_data_CVA6,
                               estimate_omega = 1L,
                               omega_prior_sd = sero_omega_prior_sd,
                               sigma_by_state = 0L)),
  # SIS
  tar_stan_mcmc(
    name = "CVA6_sr_shared",
    stan_files = "Stan/antibody_sero_sis.stan",
    data = sero_sr_shared_data_CVA6,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    cpp_options = stan_cpp, threads_per_chain = sero_threads,
    variables = sero_vars$sis,
    seed = sero_seed, refresh = sero_refresh
  ),
  
  # SIIS
  tar_stan_mcmc(
    name = "CVA6_sr_3_shared",
    stan_files = "Stan/antibody_sero_siis.stan",
    data = sero_sr_shared_data_CVA6,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    cpp_options = stan_cpp, threads_per_chain = sero_threads,
    variables = sero_vars$siis,
    seed = sero_seed, refresh = sero_refresh
  ),
  
  # EV68
  tar_target(sero_stan_data_EV68_shared,
             make_sero_stan_data(ab_prep_EV68$serum_meta, log_phi_draws_EV68,
                                 sigma_by_state = 0L)),
  # SI
  tar_stan_mcmc(
    name = "EV68_shared",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_stan_data_EV68_shared,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$si,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
  ),
  
  # SII
  tar_stan_mcmc(
    name = "EV68_three_comp_shared",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_stan_data_EV68_shared,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$sii,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_target(sero_sr_shared_data_EV68,
             add_sero_switches(sero_stan_data_EV68,
                               estimate_omega = 1L,
                               omega_prior_sd = sero_omega_prior_sd,
                               sigma_by_state = 0L)),
  # SIS
  tar_stan_mcmc(
    name = "EV68_sr_shared",
    stan_files = "Stan/antibody_sero_sis.stan",
    data = sero_sr_shared_data_EV68,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    cpp_options = stan_cpp, threads_per_chain = sero_threads,
    variables = sero_vars$sis, seed = sero_seed, refresh = sero_refresh
  ),
  
  # SIIS
  tar_stan_mcmc(
    name = "EV68_sr_3_shared",
    stan_files = "Stan/antibody_sero_siis.stan",
    data = sero_sr_shared_data_EV68,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    cpp_options = stan_cpp, threads_per_chain = sero_threads,
    variables = sero_vars$siis, seed = sero_seed, refresh = sero_refresh
  ),
  
  
  # EV71
  tar_target(sero_stan_data_EV71_shared,
             make_sero_stan_data(ab_prep_EV71$serum_meta, log_phi_draws_EV71,
                                 sigma_by_state = 0L)),
  
  # SII
  tar_stan_mcmc(
    name = "EV71_shared",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_stan_data_EV71_shared,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$sii,
    cpp_options = stan_cpp,
    threads_per_chain = sero_threads
  ),
  
  # SI
  tar_stan_mcmc(
    name = "EV71_two_comp_shared",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_stan_data_EV71_shared,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$si,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_target(sero_sr_shared_data_EV71,
             add_sero_switches(sero_stan_data_EV71,
                               estimate_omega = 1L,
                               omega_prior_sd = sero_omega_prior_sd,
                               sigma_by_state = 0L)),
  
  # SIS
  tar_stan_mcmc(
    name = "EV71_sr_shared",
    stan_files = "Stan/antibody_sero_sis.stan",
    data = sero_sr_shared_data_EV71,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    cpp_options = stan_cpp, threads_per_chain = sero_threads,
    variables = sero_vars$sis, seed = sero_seed, refresh = sero_refresh
  ),
  
  # SIIS
  tar_stan_mcmc(
    name = "EV71_sr_3_shared",
    stan_files = "Stan/antibody_sero_siis.stan",
    data = sero_sr_shared_data_EV71,
    iter_warmup = sero_sr_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    adapt_delta = sero_sr_adapt_delta,
    cpp_options = stan_cpp, threads_per_chain = sero_threads,
    variables = sero_vars$siis, seed = sero_seed, refresh = sero_refresh
  ),
  
  
  tar_target(
    model_grid_loo,
    dplyr::bind_rows(
      compare_sero_models_loo(list(
        SI_shared    = CVA6_two_comp_shared_draws_antibody_sero_two_compartment,
        SI_bystate   = CVA6_two_comp_draws_antibody_sero_two_compartment,
        SII_shared   = CVA6_shared_draws_antibody_sero,
        SII_bystate  = CVA6_draws_antibody_sero,
        SIS_shared   = CVA6_sr_shared_draws_antibody_sero_sis,
        SIIS_shared  = CVA6_sr_3_shared_draws_antibody_sero_siis
      ), virus = "CVA6"),
      compare_sero_models_loo(list(
        SI_shared    = EV71_two_comp_shared_draws_antibody_sero_two_compartment,
        SI_bystate   = EV71_two_comp_draws_antibody_sero_two_compartment,
        SII_shared   = EV71_shared_draws_antibody_sero,
        SII_bystate  = EV71_draws_antibody_sero,
        SIS_shared   = EV71_sr_shared_draws_antibody_sero_sis,
        SIIS_shared  = EV71_sr_3_shared_draws_antibody_sero_siis
      ), virus = "EV71"),
      compare_sero_models_loo(list(
        SI_shared    = EV68_shared_draws_antibody_sero_two_compartment,
        SI_bystate   = EV68_draws_antibody_sero_two_compartment,
        SII_shared   = EV68_three_comp_shared_draws_antibody_sero,
        SIS_shared   = EV68_sr_shared_draws_antibody_sero_sis
      ), virus = "EV68")
    )
  ),
  
  #### convergence diagnostics ####
  tar_target(
    convergence_df,
    label_convergence_failures(dplyr::bind_rows(
      # --- main fits, by-state sigma ---
      sero_convergence_row(CVA6_summary_antibody_sero,
                           CVA6_draws_antibody_sero, "CVA6", "SII"),
      sero_convergence_row(EV71_summary_antibody_sero,
                           EV71_draws_antibody_sero, "EV71", "SII"),
      sero_convergence_row(EV68_three_comp_summary_antibody_sero,
                           EV68_three_comp_draws_antibody_sero, "EV68", "SII"),
      sero_convergence_row(CVA6_two_comp_summary_antibody_sero_two_compartment,
                           CVA6_two_comp_draws_antibody_sero_two_compartment, "CVA6", "SI"),
      sero_convergence_row(EV71_two_comp_summary_antibody_sero_two_compartment,
                           EV71_two_comp_draws_antibody_sero_two_compartment, "EV71", "SI"),
      sero_convergence_row(EV68_summary_antibody_sero_two_compartment,
                           EV68_draws_antibody_sero_two_compartment, "EV68", "SI"),
      
      # --- shared sigma ---
      sero_convergence_row(CVA6_shared_summary_antibody_sero,
                           CVA6_shared_draws_antibody_sero, "CVA6", "SII", "shared"),
      sero_convergence_row(EV71_shared_summary_antibody_sero,
                           EV71_shared_draws_antibody_sero, "EV71", "SII", "shared"),
      sero_convergence_row(EV68_three_comp_shared_summary_antibody_sero,
                           EV68_three_comp_shared_draws_antibody_sero, "EV68", "SII", "shared"),
      sero_convergence_row(CVA6_two_comp_shared_summary_antibody_sero_two_compartment,
                           CVA6_two_comp_shared_draws_antibody_sero_two_compartment,
                           "CVA6", "SI", "shared"),
      sero_convergence_row(EV71_two_comp_shared_summary_antibody_sero_two_compartment,
                           EV71_two_comp_shared_draws_antibody_sero_two_compartment,
                           "EV71", "SI", "shared"),
      sero_convergence_row(EV68_shared_summary_antibody_sero_two_compartment,
                           EV68_shared_draws_antibody_sero_two_compartment, "EV68", "SI", "shared"),
      
      # --- seroreversion, by-state sigma ---
      sero_convergence_row(CVA6_sr_summary_antibody_sero_sis,
                           CVA6_sr_draws_antibody_sero_sis, "CVA6", "SIS"),
      sero_convergence_row(EV71_sr_summary_antibody_sero_sis,
                           EV71_sr_draws_antibody_sero_sis, "EV71", "SIS"),
      sero_convergence_row(EV68_sr_summary_antibody_sero_sis,
                           EV68_sr_draws_antibody_sero_sis, "EV68", "SIS"),
      sero_convergence_row(CVA6_sr_3_summary_antibody_sero_siis,
                           CVA6_sr_3_draws_antibody_sero_siis, "CVA6", "SIIS"),
      sero_convergence_row(EV71_sr_3_summary_antibody_sero_siis,
                           EV71_sr_3_draws_antibody_sero_siis, "EV71", "SIIS"),
      sero_convergence_row(EV68_sr_3_summary_antibody_sero_siis,
                           EV68_sr_3_draws_antibody_sero_siis, "EV68", "SIIS"),
      
      # --- seroreversion, shared sigma (isolate omega from sigma) ---
      sero_convergence_row(CVA6_sr_shared_summary_antibody_sero_sis,
                           CVA6_sr_shared_draws_antibody_sero_sis, "CVA6", "SIS", "shared"),
      sero_convergence_row(EV71_sr_shared_summary_antibody_sero_sis,
                           EV71_sr_shared_draws_antibody_sero_sis, "EV71", "SIS", "shared"),
      sero_convergence_row(EV68_sr_shared_summary_antibody_sero_sis,
                           EV68_sr_shared_draws_antibody_sero_sis, "EV68", "SIS", "shared"),
      sero_convergence_row(CVA6_sr_3_shared_summary_antibody_sero_siis,
                           CVA6_sr_3_shared_draws_antibody_sero_siis, "CVA6", "SIIS", "shared"),
      sero_convergence_row(EV71_sr_3_shared_summary_antibody_sero_siis,
                           EV71_sr_3_shared_draws_antibody_sero_siis, "EV71", "SIIS", "shared"),
      sero_convergence_row(EV68_sr_3_shared_summary_antibody_sero_siis,
                           EV68_sr_3_shared_draws_antibody_sero_siis, "EV68", "SIIS", "shared")
    ))
  ),
  
  tar_target(
    file_convergence_csv,
    {
      ensure_dir("outputs/diagnostics")
      write.csv(convergence_df, "outputs/diagnostics/convergence_summary.csv",
                row.names = FALSE)
      "outputs/diagnostics/convergence_summary.csv"
    },
    format = "file"
  ),
  
  tar_target(
    file_convergence_table,
    save_gt_pdf(make_convergence_table(convergence_df),
                "outputs/diagnostics/convergence_summary.pdf"),
    format = "file"
  ),
  
  #### stratification by survey year ####
  
  ## CVA6
  tar_target(sero_data_CVA6_2006,
             subset_sero_data_by_year(sero_stan_data_CVA6_shared,
                                      ab_prep_CVA6$serum_meta, 2006)),
  tar_target(sero_data_CVA6_2011,
             subset_sero_data_by_year(sero_stan_data_CVA6_shared,
                                      ab_prep_CVA6$serum_meta, 2011)),
  tar_target(sero_data_CVA6_2017,
             subset_sero_data_by_year(sero_stan_data_CVA6_shared,
                                      ab_prep_CVA6$serum_meta, 2017)),
  
  tar_stan_mcmc(
    name = "CVA6_yr2006",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_data_CVA6_2006,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$sii,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_stan_mcmc(
    name = "CVA6_yr2011",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_data_CVA6_2011,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$sii,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_stan_mcmc(
    name = "CVA6_yr2017",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_data_CVA6_2017,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$sii,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_target(
    CVA6_year_contrasts,
    dplyr::mutate(
      strata_contrast_all(list(
        `2006` = CVA6_yr2006_draws_antibody_sero,
        `2011` = CVA6_yr2011_draws_antibody_sero,
        `2017` = CVA6_yr2017_draws_antibody_sero
      )),
      virus = "CVA6")
  ),
  
  tar_target(
    CVA6_year_trend,
    dplyr::bind_rows(lapply(c("lambda_1", "kappa"), function(p)
      strata_trend(list(
        CVA6_yr2006_draws_antibody_sero,
        CVA6_yr2011_draws_antibody_sero,
        CVA6_yr2017_draws_antibody_sero
      ), par = p)))
  ),
  
  tar_target(
    CVA6_year_seroprev,
    dplyr::bind_rows(
      strata_seroprev(CVA6_yr2006_draws_antibody_sero,
                      sero_data_CVA6_2006$age_max, "2006"),
      strata_seroprev(CVA6_yr2011_draws_antibody_sero,
                      sero_data_CVA6_2011$age_max, "2011"),
      strata_seroprev(CVA6_yr2017_draws_antibody_sero,
                      sero_data_CVA6_2017$age_max, "2017"),
      strata_seroprev(CVA6_shared_draws_antibody_sero,
                      sero_plot_data_CVA6$age_max, "pooled"))
  ),
  
  tar_target(
    file_CVA6_by_year,
    {
      ensure_dir("outputs/stratification")
      write.csv(CVA6_year_contrasts,
                "outputs/stratification/CVA6_by_year_contrasts.csv", row.names = FALSE)
      write.csv(CVA6_year_trend,
                "outputs/stratification/CVA6_by_year_trend.csv", row.names = FALSE)
      write.csv(CVA6_year_seroprev,
                "outputs/stratification/CVA6_by_year_seroprevalence.csv", row.names = FALSE)
      c("outputs/stratification/CVA6_by_year_contrasts.csv",
        "outputs/stratification/CVA6_by_year_trend.csv",
        "outputs/stratification/CVA6_by_year_seroprevalence.csv")
    },
    format = "file"
  ),
  
  tar_target(
    file_CVA6_by_year_foi,
    {
      ensure_dir("outputs/stratification")
      g <- plot_strata_foi(
        virus="CVA6",
        
        list(
          `2006`  = CVA6_yr2006_draws_antibody_sero,
          `2011`  = CVA6_yr2011_draws_antibody_sero,
          `2017`  = CVA6_yr2017_draws_antibody_sero,
          pooled  = CVA6_shared_draws_antibody_sero
        ), 
        
        sero_plot_data_CVA6$age_max)
      
      save_plot_pdf(g, "outputs/stratification/CVA6_by_year_FOI.pdf", 8, 6)
    },
    format = "file"
  ),
  
  
  ## EV71
  tar_target(sero_data_EV71_2006,
             subset_sero_data_by_year(sero_stan_data_EV71_shared,
                                      ab_prep_EV71$serum_meta, 2006)),
  tar_target(sero_data_EV71_2011,
             subset_sero_data_by_year(sero_stan_data_EV71_shared,
                                      ab_prep_EV71$serum_meta, 2011)),
  tar_target(sero_data_EV71_2017,
             subset_sero_data_by_year(sero_stan_data_EV71_shared,
                                      ab_prep_EV71$serum_meta, 2017)),
  
  tar_stan_mcmc(
    name = "EV71_yr2006",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_data_EV71_2006,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$sii,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_stan_mcmc(
    name = "EV71_yr2011",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_data_EV71_2011,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$sii,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_stan_mcmc(
    name = "EV71_yr2017",
    stan_files = "Stan/antibody_sero.stan",
    data = sero_data_EV71_2017,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$sii,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_target(
    EV71_year_contrasts,
    dplyr::mutate(
      strata_contrast_all(list(
        `2006` = EV71_yr2006_draws_antibody_sero,
        `2011` = EV71_yr2011_draws_antibody_sero,
        `2017` = EV71_yr2017_draws_antibody_sero
      )),
      virus = "EV71")
  ),
  
  tar_target(
    EV71_year_trend,
    dplyr::bind_rows(lapply(c("lambda_1", "kappa"), function(p)
      strata_trend(list(
        EV71_yr2006_draws_antibody_sero,
        EV71_yr2011_draws_antibody_sero,
        EV71_yr2017_draws_antibody_sero
      ), par = p)))
  ),
  
  tar_target(
    EV71_year_seroprev,
    dplyr::bind_rows(
      strata_seroprev(EV71_yr2006_draws_antibody_sero,
                      sero_data_EV71_2006$age_max, "2006"),
      strata_seroprev(EV71_yr2011_draws_antibody_sero,
                      sero_data_EV71_2011$age_max, "2011"),
      strata_seroprev(EV71_yr2017_draws_antibody_sero,
                      sero_data_EV71_2017$age_max, "2017"),
      strata_seroprev(EV71_shared_draws_antibody_sero,
                      sero_plot_data_EV71$age_max, "pooled"))
  ),
  
  tar_target(
    file_EV71_by_year,
    {
      ensure_dir("outputs/stratification")
      write.csv(EV71_year_contrasts,
                "outputs/stratification/EV71_by_year_contrasts.csv", row.names = FALSE)
      write.csv(EV71_year_trend,
                "outputs/stratification/EV71_by_year_trend.csv", row.names = FALSE)
      write.csv(EV71_year_seroprev,
                "outputs/stratification/EV71_by_year_seroprevalence.csv", row.names = FALSE)
      c("outputs/stratification/EV71_by_year_contrasts.csv",
        "outputs/stratification/EV71_by_year_trend.csv",
        "outputs/stratification/EV71_by_year_seroprevalence.csv")
    },
    format = "file"
  ),
  
  tar_target(
    file_EV71_by_year_foi,
    {
      ensure_dir("outputs/stratification")
      g <- plot_strata_foi(
        virus="EV-A71",
        
        list(
        `2006`  = EV71_yr2006_draws_antibody_sero,
        `2011`  = EV71_yr2011_draws_antibody_sero,
        `2017`  = EV71_yr2017_draws_antibody_sero,
        pooled  = EV71_shared_draws_antibody_sero
        ), 
        
        sero_plot_data_EV71$age_max)
      
      save_plot_pdf(g, "outputs/stratification/EV71_by_year_FOI.pdf", 8, 6)
    },
    format = "file"
  ),
  
  ## EV68
  tar_target(sero_data_EV68_2006,
             subset_sero_data_by_year(sero_stan_data_EV68,
                                      ab_prep_EV68$serum_meta, 2006)),
  tar_target(sero_data_EV68_2011,
             subset_sero_data_by_year(sero_stan_data_EV68,
                                      ab_prep_EV68$serum_meta, 2011) %>%
               list_modify(sigma_by_state=0)), # model selection prefers shared SI for 2011 only (by-state did not converge)
  tar_target(sero_data_EV68_2017,
             subset_sero_data_by_year(sero_stan_data_EV68,
                                      ab_prep_EV68$serum_meta, 2017)),
  
  tar_stan_mcmc(
    name = "EV68_yr2006",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_data_EV68_2006,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$si,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_stan_mcmc(
    name = "EV68_yr2011",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_data_EV68_2011,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$si,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_stan_mcmc(
    name = "EV68_yr2017",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_data_EV68_2017,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$si,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_target(
    EV68_year_contrasts,
    dplyr::mutate(
      strata_contrast_all(list(
        `2006` = EV68_yr2006_draws_antibody_sero_two_compartment,
        `2011` = EV68_yr2011_draws_antibody_sero_two_compartment,
        `2017` = EV68_yr2017_draws_antibody_sero_two_compartment
      )),
      virus = "EV68")
  ),
  
  tar_target(
    EV68_year_trend,
    dplyr::bind_rows(lapply(c("lambda_1", "kappa"), function(p)
      strata_trend(list(
        EV68_yr2006_draws_antibody_sero_two_compartment,
        EV68_yr2011_draws_antibody_sero_two_compartment,
        EV68_yr2017_draws_antibody_sero_two_compartment
      ), par = p)))
  ),
  
  tar_target(
    EV68_year_seroprev,
    dplyr::bind_rows(
      strata_seroprev(EV68_yr2006_draws_antibody_sero_two_compartment,
                      sero_data_EV68_2006$age_max, "2006"),
      strata_seroprev(EV68_yr2011_draws_antibody_sero_two_compartment,
                      sero_data_EV68_2011$age_max, "2011"),
      strata_seroprev(EV68_yr2017_draws_antibody_sero_two_compartment,
                      sero_data_EV68_2017$age_max, "2017"),
      strata_seroprev(EV68_draws_antibody_sero_two_compartment,
                      sero_plot_data_EV68$age_max, "pooled"))
  ),
  
  tar_target(
    file_EV68_by_year,
    {
      ensure_dir("outputs/stratification")
      write.csv(EV68_year_contrasts,
                "outputs/stratification/EV68_by_year_contrasts.csv", row.names = FALSE)
      write.csv(EV68_year_trend,
                "outputs/stratification/EV68_by_year_trend.csv", row.names = FALSE)
      write.csv(EV68_year_seroprev,
                "outputs/stratification/EV68_by_year_seroprevalence.csv", row.names = FALSE)
      c("outputs/stratification/EV68_by_year_contrasts.csv",
        "outputs/stratification/EV68_by_year_trend.csv",
        "outputs/stratification/EV68_by_year_seroprevalence.csv")
    },
    format = "file"
  ),
  
  tar_target(
    file_EV68_by_year_foi,
    {
      ensure_dir("outputs/stratification")
      g <- plot_strata_foi(
        virus = "EV-D68",
        
        list(
        `2006`  = EV68_yr2006_draws_antibody_sero_two_compartment,
        `2011`  = EV68_yr2011_draws_antibody_sero_two_compartment,
        `2017`  = EV68_yr2017_draws_antibody_sero_two_compartment,
        pooled  = EV68_draws_antibody_sero_two_compartment
      ), 
      sero_plot_data_EV68$age_max)
      
      save_plot_pdf(g, "outputs/stratification/EV68_by_year_FOI.pdf", 8, 6)
    },
    format = "file"
  ),
  
  # pre vs post 2014 (EV-D68)
  tar_target(sero_data_EV68_pre2014,
             subset_sero_data_by_year(sero_stan_data_EV68,
                                      ab_prep_EV68$serum_meta, c(2006, 2011))),
  tar_target(sero_data_EV68_post2014,
             subset_sero_data_by_year(sero_stan_data_EV68,
                                      ab_prep_EV68$serum_meta, 2017)),
  
  tar_stan_mcmc(
    name = "EV68_pre2014",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_data_EV68_pre2014,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$si,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_stan_mcmc(
    name = "EV68_post2014",
    stan_files = "Stan/antibody_sero_two_compartment.stan",
    data = sero_data_EV68_post2014,
    iter_warmup = sero_iter_warmup,
    iter_sampling = sero_iter_sampling - sero_iter_warmup,
    chains = sero_chains, parallel_chains = sero_parallel,
    seed = sero_seed, refresh = sero_refresh,
    variables = sero_vars$si,
    cpp_options = stan_cpp, threads_per_chain = sero_threads
  ),
  
  tar_target(
    EV68_prepost_contrast,
    dplyr::mutate(
      strata_contrast(EV68_pre2014_draws_antibody_sero_two_compartment,
                      EV68_post2014_draws_antibody_sero_two_compartment,
                      "pre-2014 (2006, 2011)", "post-2014 (2017)"),
      virus = "EV68",
      n_a = sero_data_EV68_pre2014$n,
      n_b = sero_data_EV68_post2014$n)
  ),
  
  tar_target(
    EV68_prepost_seroprev,
    dplyr::bind_rows(
      strata_seroprev(EV68_pre2014_draws_antibody_sero_two_compartment,
                      sero_data_EV68_pre2014$age_max, "pre-2014"),
      strata_seroprev(EV68_post2014_draws_antibody_sero_two_compartment,
                      sero_data_EV68_post2014$age_max, "post-2014"),
      strata_seroprev(EV68_draws_antibody_sero_two_compartment,
                      sero_plot_data_EV68$age_max, "pooled"))
  ),
  
  tar_target(
    file_EV68_prepost,
    {
      ensure_dir("outputs/stratification")
      write.csv(EV68_prepost_contrast,
                "outputs/stratification/EV68_pre_post_2014_contrast.csv",
                row.names = FALSE)
      write.csv(EV68_prepost_seroprev,
                "outputs/stratification/EV68_pre_post_2014_seroprevalence.csv",
                row.names = FALSE)
      c("outputs/stratification/EV68_pre_post_2014_contrast.csv",
        "outputs/stratification/EV68_pre_post_2014_seroprevalence.csv")
    },
    format = "file"
  ),
  
  tar_target(
    file_EV68_prepost_foi,
    {
      ensure_dir("outputs/stratification")
      g <- plot_strata_foi(
        virus = "EV-D68",
        
        list(
        `pre-2014`  = EV68_pre2014_draws_antibody_sero_two_compartment,
        `post-2014` = EV68_post2014_draws_antibody_sero_two_compartment,
        pooled      = EV68_draws_antibody_sero_two_compartment
      ), 
      
      sero_plot_data_EV68$age_max)
      
      save_plot_pdf(g, "outputs/stratification/EV68_pre_post_2014_FOI.pdf", 8, 6)
    },
    format = "file"
  )
)
