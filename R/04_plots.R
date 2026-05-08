library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gt)
library(matrixStats)
library(viridis)

# Global palettes / labels

serostate_cols <- c(
  S   = "#1F77B4FF",   # blue
  Ipp = "#FF7F0EFF",   # orange
  Ip  = "#8B4513"      # brown
)

serostate_levels <- c("S", "Ipp", "Ip")
serostate_labels <- c(S = "S", Ipp = "I++", Ip = "I+")

virus_cols <- c(
  EV71 = "#A020F0",
  CVA6 = "#009E73",
  EV68 = "#D4AF37"
)
virus_labels <- c(
  EV71 = "EV-A71",
  CVA6 = "CVA6",
  EV68 = "EV-D68"
)

scale_virus_color <- function(...) {
  scale_color_manual(values = virus_cols, labels = virus_labels, ...)
}

scale_virus_fill <- function(...) {
  scale_fill_manual(values = virus_cols, labels = virus_labels, ...)
}

scale_serostate_color <- function(state_names = serostate_levels, ...) {
  scale_color_manual(
    values = serostate_cols[state_names],
    breaks = state_names,
    labels = serostate_labels[state_names],
    ...
  )
}

scale_serostate_fill <- function(state_names = serostate_levels, ...) {
  scale_fill_manual(
    values = serostate_cols[state_names],
    breaks = state_names,
    labels = serostate_labels[state_names],
    ...
  )
}

# distinguish between 3-state model (EV71, CVA6) and 2-state model
detect_n_states <- function(summary_df) {
  if (any(grepl("^mu\\[3\\]", summary_df$variable))) {
    3
  } else {
    2
  }
}

# extract prob_by_group mean matrix from cmdstan summary
extract_prob_matrix <- function(summary_df, n, k = 3) {
  probs <- summary_df %>%
    filter(str_detect(variable, "^prob_by_group\\[")) %>%
    select(variable, mean)
  
  ij <- stringr::str_match(probs$variable, "^prob_by_group\\[(\\d+),(\\d+)\\]$")
  i <- as.integer(ij[, 2])
  j <- as.integer(ij[, 3])
  
  mat <- matrix(NA_real_, nrow = n, ncol = k)
  mat[cbind(i, j)] <- probs$mean
  mat
}

# FOI plot: "short" from (lambda_1,beta) draws, "long" from lambda_long summary
plot_foi_short_long <- function(virus, plot_data, summary_df, draws_df) {
  max_age <- plot_data$age_max
  age_seq <- 0:(max_age - 1)
  
  post <- draws_df %>% select(lambda_1, beta)
  
  FOI_mat <- vapply(
    seq_len(nrow(post)),
    function(r) post$lambda_1[r] * exp(-post$beta[r] * age_seq),
    numeric(length(age_seq))
  )
  FOI_mat <- t(FOI_mat) # [draw x age]
  
  short_df <- data.frame(
    Age = age_seq,
    mean = colMeans(FOI_mat),
    q5 = apply(FOI_mat, 2, quantile, 0.05),
    q95 = apply(FOI_mat, 2, quantile, 0.95)
  ) %>% filter(Age > 0)
  
  long_df <- summary_df %>%
    filter(str_detect(variable, "^lambda_long\\[")) %>%
    transmute(
      Age = as.integer(str_match(variable, "^lambda_long\\[(\\d+)\\]$")[, 2]) - 1L,
      mean = mean, q5 = q5, q95 = q95
    ) %>% filter(Age > 0)
  
  # short_col <- if (virus == "EV71") "#A020F0" else "#009E73"
  short_col <- virus_cols[[virus]]
  
  ggplot() +
    geom_ribbon(data = long_df, aes(x = Age, ymin = q5, ymax = q95, fill = "Long"), alpha = 0.2) +
    geom_line(data = long_df, aes(x = Age, y = mean, color = "Long"), linewidth = 1.2) +
    geom_ribbon(data = short_df, aes(x = Age, ymin = q5, ymax = q95, fill = "Short"), alpha = 0.2) +
    geom_line(data = short_df, aes(x = Age, y = mean, color = "Short"), linewidth = 1.2) +
    scale_color_manual(values = c(Long = "black", Short = short_col)) +
    scale_fill_manual(values = c(Long = "black", Short = short_col)) +
    labs(
      title = paste0("FoI Posterior Estimates: ", pretty_virus(virus)),
      x = "Age (Years)", y = "Force of Infection", color = NULL, fill = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_x_log10()
}

plot_foi_combined <- function(draws_combined_melt, max_age = 92) {
  viruses <- c("EV71", "CVA6", "EV68")
  age_seq <- 1:max_age
  
  dm <- draws_combined_melt %>%
    filter(Virus %in% viruses, variable %in% c("lambda_1", "beta")) %>%
    group_by(Virus, variable) %>%
    mutate(draw = row_number()) %>%
    ungroup()
  
  wide <- dm %>%
    pivot_wider(
      id_cols = setdiff(names(dm), c("variable", "value")),
      names_from = variable,
      values_from = value
    ) %>%
    # defensively coerce away list-cols if they exist
    mutate(across(c(lambda_1, beta), ~ as.numeric(unlist(.x))))
  
  foi_summ <- wide %>%
    tidyr::crossing(Age = age_seq) %>%
    mutate(FOI = lambda_1 * exp(-beta * Age)) %>%
    group_by(Virus, Age) %>%
    summarise(
      mean = mean(FOI, na.rm = TRUE),
      q5 = quantile(FOI, 0.05, na.rm = TRUE),
      q95 = quantile(FOI, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(foi_summ, aes(x = Age, y = mean, color = Virus, fill = Virus, group = Virus)) +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    scale_virus_color() +
    scale_virus_fill() +
    labs(x = "Age (Years)", y = "Force of Infection") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_x_log10()
}

plot_param_densities_combined <- function(draws_combined_melt) {
  keep <- c("mu[1]", "mu[2]", "mu[3]", "sigma", "psi", "lambda_1", "beta")
  dd <- draws_combined_melt %>% filter(variable %in% keep) %>% filter(!(variable == "beta" & value>0.8))
  
  var_labs <- c(
    "mu[1]" = "mu[S]",
    "mu[2]" = "mu[I[\"++\"]]",
    "mu[3]" = "mu[I[\"+\"]]",
    "sigma" = "sigma",
    "psi" = "psi",
    "lambda_1" = "lambda[1]",
    "beta" = "kappa"
  )
  
  dd$Virus <- factor(dd$Virus, levels = c("CVA6", "EV71", "EV68"))
  
  ggplot(dd, aes(x = value, color = Virus, fill = Virus)) +
    geom_density(alpha = 0.5) +
    facet_wrap(
      ~variable,
      scales = "free",
      labeller = as_labeller(var_labs, label_parsed)
    ) +
    scale_virus_color() +
    scale_virus_fill() +
    labs(x = "Value", y = "Density") +
    theme_minimal() +
    theme(
      legend.position = c(0.75, 0.025),
      legend.justification = c("right", "bottom"),
      strip.text = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      legend.key.height = grid::unit(0.75, "cm"),
      legend.key.width = grid::unit(0.75, "cm"),
      legend.key.spacing.y = grid::unit(0.25, "cm")
    )
}

make_param_table <- function(virus, summary_df) {
  n_states <- detect_n_states(summary_df)
  
  if (n_states == 3) {
    vars <- c("mu[1]", "mu[3]", "mu[2]", "sigma", "psi", "lambda_1", "beta")
    recode_map <- c(
      "mu[1]" = "\\mu_{S}",
      "mu[3]" = "\\mu_{I_{++}}",
      "mu[2]" = "\\mu_{I_{+}}",
      "sigma" = "\\sigma",
      "psi" = "\\psi",
      "lambda_1" = "\\lambda_{1}",
      "beta" = "\\kappa"
    )
    desired_order <- c(
      "\\mu_{S}",
      "\\mu_{I_{++}}",
      "\\mu_{I_{+}}",
      "\\sigma",
      "\\psi",
      "\\lambda_{1}",
      "\\kappa"
    )
  } else {
    vars <- c("mu[1]", "mu[2]", "sigma", "lambda_1", "beta")
    recode_map <- c(
      "mu[1]" = "\\mu_{S}",
      "mu[2]" = "\\mu_{I_{++}}",
      "sigma" = "\\sigma",
      "lambda_1" = "\\lambda_{1}",
      "beta" = "\\kappa"
    )
    desired_order <- c(
      "\\mu_{S}",
      "\\mu_{I_{++}}",
      "\\sigma",
      "\\lambda_{1}",
      "\\kappa"
    )
  }
  
  
  tab <- summary_df %>%
    filter(variable %in% vars) %>%
    select(variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail) %>%
    mutate(variable = recode(variable, !!!recode_map))
  
  tab$variable <- factor(tab$variable, levels = desired_order)
  tab <- tab %>% arrange(variable)
  tab$variable <- paste0("$", tab$variable, "$")
  
  gt(tab) %>%
    tab_header(
      title = paste0(pretty_virus(virus), ": Summary Statistics of Fitted Parameters")
    ) %>%
    gt_theme_nytimes() %>%
    fmt_markdown(columns = variable) %>%
    fmt_number(columns = c(-variable), decimals = 3) %>%
    fmt_number(columns = c(ess_bulk, ess_tail), decimals = 0) %>%
    cols_align(
      align = "center",
      columns = variable
    ) %>%
    cols_align(
      align = "right",
      columns = c(mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)
    ) %>%
    tab_style(
      style = list(
        cell_text(weight = "bold"),
        cell_borders(sides = "bottom", color = "#222222", weight = px(1))
      ),
      locations = cells_column_labels(everything())
    ) %>%
    cols_width(
      variable ~ px(105),
      mean ~ px(68), median ~ px(68), sd ~ px(68),
      q5 ~ px(68), q95 ~ px(68), rhat ~ px(58),
      ess_bulk ~ px(72), ess_tail ~ px(72)
    ) %>%
    tab_options(
      table.font.size = px(14),          
      heading.title.font.size = px(14),
      heading.align = "center",
      column_labels.font.weight = "bold",
      column_labels.padding = px(5),
      data_row.padding = px(4),
      table.border.top.color = "#222222",
      table.border.top.width = px(1),
      table.border.bottom.color = "#222222",
      table.border.bottom.width = px(1),
      column_labels.border.bottom.color = "#222222",
      column_labels.border.bottom.width = px(1),
      table_body.hlines.color = "#D9D9D9",
      table_body.hlines.width = px(0.5)
    )
}


plot_phi_points <- function(virus, plot_data, summary_df, eps = 0.95) {
  n <- plot_data$n
  ages <- plot_data$ages
  phi <- plot_data$phi_median
  
  prob_mat <- extract_prob_matrix(summary_df, n = n, k = 3)
  
  # sample a state for visualization (consistent with your earlier approach)
  set.seed(1)
  state <- integer(n)
  below <- logical(n)
  for (i in seq_len(n)) {
    p <- prob_mat[i, ]
    p[is.na(p)] <- 0
    p <- p / sum(p)
    s <- sample.int(3, 1, prob = p)
    state[i] <- s
    below[i] <- (p[s] < eps)
  }
  
  df <- data.frame(
    age = ages,
    phi = phi,
    state = factor(state, levels = 1:3, labels = serostate_levels),
    conf = ifelse(below, "<95%", ">95%")
  ) %>% mutate(key = interaction(state, conf, sep = "_"))
  
  ggplot(df, aes(x = age, y = phi, color = state, shape = conf)) +
    geom_point(alpha = 0.9) +
    scale_y_log10() +
    scale_serostate_color(state_names = serostate_levels) +
    scale_shape_manual(values = c(">95%" = 16, "<95%" = 0)) +
    labs(
      title = paste0("Antibody concentration vs Age: ", pretty_virus(virus)),
      x = "Age (Years)", y = expression(phi),
      color = "State", shape = "Confidence"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.title.y = element_text(
        face = "bold",
        size = 16,
        angle = 0,
        vjust = 0.5
      )
    )
}

simulate_serodynamics <- function(lambda_long, psi, age_max) {
  # returns matrix [age_max x 3] columns S, I, R for ages 0..age_max-1
  S <- numeric(age_max); I <- numeric(age_max); R <- numeric(age_max)
  S[1] <- 1; I[1] <- 0; R[1] <- 0
  
  for (a in 2:age_max) {
    lam <- lambda_long[a]
    oldS <- S[a - 1]
    oldI <- I[a - 1]
    
    S[a] <- oldS * exp(-lam)
    
    if (abs(psi - lam) > 1e-12) {
      I[a] <- exp(-psi) * oldI + (lam / (psi - lam)) * oldS * (exp(-lam) - exp(-psi))
    } else {
      I[a] <- exp(-psi) * (oldI + psi * oldS)
    }
    
    R[a] <- 1 - (S[a] + I[a])
  }
  
  cbind(S = S, I = I, R = R)
}

plot_weighted_avg_phi <- function(
    virus,
    plot_data,
    draws_df,
    summary_df,
    n_samples = 1000,
    seed = 1,
    eps = 0.95
) {
  
  n_states <- detect_n_states(summary_df)
  
  age_max <- plot_data$age_max
  ages <- plot_data$ages
  tps <- 0:(age_max - 1)
  
  phi_points <- plot_data$phi_median
  n <- length(phi_points)
  
  prob_mat <- extract_prob_matrix(summary_df, n = n, k = n_states)
  
  set.seed(seed)
  state <- integer(n)
  below <- logical(n)
  
  for (i in seq_len(n)) {
    p <- prob_mat[i, ]
    p[is.na(p)] <- 0
    p <- p / sum(p)
    s <- sample.int(n_states, 1, prob = p)
    state[i] <- s
    below[i] <- (p[s] < eps)
  }
  
  state_labels <- serostate_levels
  
  pts <- data.frame(
    age = ages,
    phi = phi_points,
    state = factor(state_labels[state], levels = state_labels),
    conf = ifelse(below, "<95%", ">95%")
  )
  
  # Compute state-specific mean φ from (mu_x, sigma)
  # E[lognormal] = exp(mu + sigma^2 / 2)
  
  mu_x <- summary_df %>%
    filter(grepl("^mu_x\\[", variable)) %>%
    arrange(variable) %>%
    pull(mean)
  
  sigma <- summary_df %>%
    filter(variable == "sigma") %>%
    pull(mean)
  
  phi_state <- exp(mu_x + 0.5 * sigma^2)
  
  if (n_states == 2) {
    phi_state <- c(phi_state, 0)
  }
  
  set.seed(seed)
  idx <- sample.int(nrow(draws_df), size = min(n_samples, nrow(draws_df)))
  
  lam_cols <- paste0("lambda_long[", 1:age_max, "]")
  stopifnot(all(lam_cols %in% names(draws_df)))
  
  sims <- lapply(idx, function(r) {
    lambda_long <- as.numeric(draws_df[r, lam_cols])
    
    if ("psi" %in% names(draws_df)) {
      simulate_serodynamics(lambda_long, draws_df$psi[r], age_max)
    } else {
      simulate_serodynamics(lambda_long, psi = 0, age_max)
    }
  })
  
  S_mat <- do.call(cbind, lapply(sims, function(M) M[, "S"]))
  I_mat <- do.call(cbind, lapply(sims, function(M) M[, "I"]))
  
  if (n_states == 3) {
    R_mat <- do.call(cbind, lapply(sims, function(M) M[, "R"]))
    probs_age <- cbind(
      S = rowMeans(S_mat),
      Ipp = rowMeans(I_mat),
      Ip = rowMeans(R_mat)
    )
  } else {
    probs_age <- cbind(
      S = rowMeans(S_mat),
      Ipp = rowMeans(I_mat),
      Ip = 0
    )
  }
  
  weighted_phi <- rowSums(
    probs_age * matrix(phi_state, nrow = nrow(probs_age), ncol = 3, byrow = TRUE)
  )
  
  line_df <- data.frame(
    Age = tps,
    weighted_phi = weighted_phi
  )
  
  ggplot() +
    geom_point(
      data = pts,
      aes(x = age, y = phi, color = state, shape = conf),
      alpha = 0.85,
      show.legend=TRUE
    ) +
    geom_line(
      data = line_df,
      aes(x = Age, y = weighted_phi),
      linewidth = 1.4
    ) +
    scale_y_log10(
      breaks = c(10,30,100,300,1000,3000)
    ) +
    scale_serostate_color(state_names = state_labels, drop = FALSE) +
    scale_shape_manual(values = c(">95%" = 16, "<95%" = 0)) +
    scale_x_continuous(
      breaks = seq(0, 90, by=10),
      limits = c(0, age_max)
    ) +
    labs(
      title = pretty_virus(virus),
      x = "Age (Years)",
      y = expression(phi),
      color = "Serostate"
      ,shape = "Confidence"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(
        size = 36,
        angle = 0,
        vjust = 0.5,
        face = "bold"
      ),
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14)
    ) + 
    guides(shape = "none",
           color = guide_legend(override.aes = list(size=7)))
}


plot_serodynamics <- function(
    virus,
    plot_data,
    draws_df,
    n_samples = 300,
    seed = 1
) {
  
  age_max <- plot_data$age_max
  ages <- 0:(age_max - 1)
  
  set.seed(seed)
  idx <- sample.int(nrow(draws_df), size = min(n_samples, nrow(draws_df)))
  
  lam_cols <- paste0("lambda_long[", 1:age_max, "]")
  stopifnot(all(lam_cols %in% names(draws_df)))
  
  # Simulate trajectories
  sims <- lapply(idx, function(r) {
    lambda_long <- as.numeric(draws_df[r, lam_cols])
    
    if ("psi" %in% names(draws_df)) {
      simulate_serodynamics(lambda_long, draws_df$psi[r], age_max)
    } else {
      simulate_serodynamics(lambda_long, psi = 0, age_max)
    }
  })
  
  S_mat <- do.call(cbind, lapply(sims, function(M) M[, "S"]))
  I_mat <- do.call(cbind, lapply(sims, function(M) M[, "I"]))
  
  if (ncol(sims[[1]]) == 3) {
    R_mat <- do.call(cbind, lapply(sims, function(M) M[, "R"]))
    state_names <- c("S", "Ipp", "Ip")
  } else {
    state_names <- c("S", "Ipp")
  }
  
  # Build summary dataframe
  df <- data.frame(Age = ages)
  
  df$S_mean <- rowMeans(S_mat)
  df$S_lo <- apply(S_mat, 1, quantile, 0.025)
  df$S_hi <- apply(S_mat, 1, quantile, 0.975)
  
  df$Ipp_mean <- rowMeans(I_mat)
  df$Ipp_lo <- apply(I_mat, 1, quantile, 0.025)
  df$Ipp_hi <- apply(I_mat, 1, quantile, 0.975)
  
  if ("Ip" %in% state_names) {
    df$Ip_mean <- rowMeans(R_mat)
    df$Ip_lo <- apply(R_mat, 1, quantile, 0.025)
    df$Ip_hi <- apply(R_mat, 1, quantile, 0.975)
  }
  
  # Convert to long format (clean mapping)
  df_long <- tidyr::pivot_longer(
    df,
    cols = -Age,
    names_to = c("State", ".value"),
    names_pattern = "(.*)_(mean|lo|hi)"
  )
  
  df_long$State <- factor(df_long$State, levels = state_names)
  
  
  # filter out (unmodelled) I+ state for EV68 (but keep legend)
  if (virus == "EV68") {
    df_long <- df_long %>% dplyr::filter(State != "Ip")
  }
  
  # make dummy df to keep consistent legend across virus plots
  legend_df <- data.frame(
    Age = NA_real_,
    mean = NA_real_,
    lo = NA_real_,
    hi = NA_real_,
    State = factor(state_names, levels = state_names)
  )
  
  ggplot(df_long, aes(x = Age)) +
    geom_ribbon(
      aes(ymin = lo, ymax = hi, fill = State),
      alpha = 0.18
    ) +
    geom_line(
      aes(y = mean, color = State),
      linewidth = 1.3
    ) +
    
    # for legend only
    geom_ribbon(
      data = legend_df,
      aes(ymin = lo, ymax = hi, fill = State),
      alpha = 0.18,
      show.legend = TRUE,
      na.rm = TRUE
    ) +
    
    # for legend only
    geom_line(
      data = legend_df,
      aes(y = mean, color = State),
      linewidth = 1.3,
      show.legend = TRUE,
      na.rm = TRUE
    ) +
    
    scale_serostate_color(state_names = state_names,
                          limits = state_names,
                          drop = FALSE) +
    scale_serostate_fill(state_names = state_names,
                         limits = state_names,
                         drop = FALSE) +
    scale_x_continuous(
      breaks = seq(0, 90, by = 10),
      limits = c(0, age_max)
    ) +
    labs(
      title = pretty_virus(virus),
      x = "Age (Years)",
      y = "Proportion",
      color = "Serostate",
      fill = "Serostate"
    ) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(
            size = 16,
            face = "bold",
            hjust = 0.5
            ),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.box = "vertical",
          legend.direction = "vertical",
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14),
          legend.key.height = grid::unit(1, "cm"),
          legend.key.width = grid::unit(2, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))
}



# new vs old FOI comparison
plot_foi_new_old <- function(
    virus,
    plot_data,
    summary_df,
    draws_df
) {
  max_age <- plot_data$age_max
  age_seq <- 0:(max_age - 1)
  
  post <- draws_df %>% dplyr::select(lambda_1, beta)
  
  FOI_mat <- vapply(
    seq_len(nrow(post)),
    function(r) post$lambda_1[r] * exp(-post$beta[r] * age_seq),
    numeric(length(age_seq))
  )
  FOI_mat <- t(FOI_mat)
  
  new_df <- data.frame(
    Age = age_seq,
    mean = colMeans(FOI_mat),
    q5 = apply(FOI_mat, 2, quantile, 0.05),
    q95 = apply(FOI_mat, 2, quantile, 0.95)
  ) %>% dplyr::filter(Age > 0)
  
  old_fit <- readRDS(file.path("data", "processed", paste0(virus, "_model5_fit.rds")))
  
  old_post <- rstan::extract(old_fit, pars = c("lambda", "beta"))
  
  old_lambda <- old_post$lambda
  old_beta   <- old_post$beta
  
  old_FOI_mat <- matrix(data=NA, nrow=max_age, ncol=length(old_beta))
  
  for (i in 1:length(old_beta)) {
    old_FOI_mat[,i] <- old_lambda[i]*exp(-old_beta[i]*age_seq)
  }
  
  old_df <- data.frame(Age = age_seq,
                       mean = rowMeans(old_FOI_mat),
                       q5   = matrixStats::rowQuantiles(old_FOI_mat, probs=0.05),
                       q95  = matrixStats::rowQuantiles(old_FOI_mat, probs=0.95)
                       ) %>% dplyr::filter(Age > 0) 
  
  virus_col   <- virus_cols [[virus]]
  
  ggplot() +
    geom_ribbon(
      data = new_df,
      aes(x = Age, ymin = q5, ymax = q95, fill = "New"),
      alpha = 0.3
    ) +
    geom_line(
      data = new_df,
      aes(x = Age, y = mean, color = "New"),
      linewidth = 1.5
    ) +
    geom_ribbon(
      data = old_df,
      aes(x = Age, ymin = q5, ymax = q95, fill = "Old"),
      alpha = 0.3
    ) +
    geom_line(
      data = old_df,
      aes(x = Age, y = mean, color = "Old"),
      linewidth = 1.5
    ) +
    labs(
      title = pretty_virus(virus),
      x = "Age (Years)",
      y = "Force of Infection",
      color = "Force of Infection",
      fill = "Force of Infection"
    ) +
    scale_x_log10() +
    scale_y_continuous(
      breaks = seq(0, 0.6, by = 0.2),
      limits = c(0, 0.66)
    ) +
    scale_color_manual(
      values = c(
        "New" = virus_col,
        "Old" = "grey80"
      ),
      labels = c(
        "New" = "New",
        "Old" = "Old"
      )
    ) + 
    scale_fill_manual(
      values = c(
        "New" = virus_col,
        "Old" = "grey80"
      ),
      labels = c(
        "New" = "New",
        "Old" = "Old"
      )
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 18, face="bold"),
      legend.text = element_text(size = 16),
      strip.text.x = element_text(size=20, face="bold"),
      legend.position = "right",
      legend.key.width = unit(3, "cm")
    )  + guides(color = guide_legend(override.aes = list(linewidth=1.5)))
}


# for combined weighted avg phi, serodynamics, and new vs old FoI plots
plot_patchwork <- function(cva6_plot, ev71_plot, ev68_plot, lines=TRUE, FOI=FALSE) {
  
  if (FOI) {
    ev71_legend <- cowplot::get_legend(
      ev71_plot + theme(legend.position = "right",
                        legend.title = element_text(size = 18, face = "bold"),
                        legend.text = element_text(size = 16),
                        legend.key.height = grid::unit(1.5, "cm"),
                        legend.key.width = if (lines) {grid::unit(4, "cm")} else {grid::unit(1, "cm")}
                        )
    )
    
    patchwork::wrap_plots(
      cva6_plot + labs(x = NULL) + theme(legend.position = "none"),
      ev71_plot + labs(x = NULL, y = NULL) + theme(legend.position = "none"),
      ev68_plot + theme(legend.position = "none"),
      patchwork::wrap_elements(full = ev71_legend),
      design = "AB\nCD"
    ) &
      theme(
        legend.position = "none",
        legend.box = "vertical",
        legend.direction = "vertical"
      )
  }
  
  else {
    
    legend_plot <- cva6_plot + 
      theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.height = grid::unit(1.5, "cm"),
        legend.key.width = if (lines) {
          grid::unit(4, "cm")
        } else {
          grid::unit(1, "cm")}
      )
    
    legend_grob <- cowplot::get_legend(legend_plot)
    
    patchwork::wrap_plots(
      cva6_plot + labs(x=NULL)         + theme(legend.position = "none"),
      ev71_plot + labs(x=NULL, y=NULL) + theme(legend.position = "none"),
      ev68_plot                        + theme(legend.position = "none"),
      patchwork::wrap_elements(full = legend_grob),
      design = "AB\nCD"
    ) # &
    #   theme(
    #     legend.position = "right",
    #     legend.box = "vertical",
    #     legend.direction = "vertical",
    #     legend.title = element_text(size = 18, face = "bold"),
    #     legend.text = element_text(size = 16),
    #     legend.key.height = grid::unit(1.5, "cm"),
    #     legend.key.width = if (lines) {
    #       grid::unit(4, "cm")
    #     } else {
    #         grid::unit(1, "cm")}
    #   )
  }
}

#### Reed-Muench ####

plot_reed_muench_fit <- function(raw_df, draws_obj, seed = 1) {
  ddf <- posterior::as_draws_df(draws_obj)
  
  phi <- ddf[["phi[1]"]]
  k1  <- ddf[["k1"]]
  
  # d_grid <- exp(seq(log(min(raw_df$dilutions)), log(max(raw_df$dilutions)), length.out = 300))
  d_grid <- seq(0.75, 350, length.out=500)
  p_mat <- sapply(d_grid, function(d) stats::pnorm(k1 * (log(phi) - log(d))))
  
  fit_df <- data.frame(
    dilution = d_grid,
    lo = apply(p_mat, 2, stats::quantile, probs = 0.025),
    med = apply(p_mat, 2, stats::quantile, probs = 0.5),
    hi = apply(p_mat, 2, stats::quantile, probs = 0.975)
  )
  
  obs_df <- raw_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      survival = outcome / n_replicates,
      ci = list(stats::binom.test(outcome, n_replicates)$conf.int),
      lo = ci[[1]],
      hi = ci[[2]]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-ci)
  
  rm_endpoint <- reed_muench_endpoint(raw_df)
  bayes_endpoint <- stats::median(phi)
  
  ggplot() +
    geom_ribbon(
      data = fit_df,
      aes(x = dilution, ymin = lo, ymax = hi),
      fill = "steelblue",
      alpha = 0.25
    ) +
    geom_line(
      data = fit_df,
      aes(x = dilution, y = med),
      linetype = "dashed",
      linewidth = 0.9
    ) +
    geom_errorbar(
      data = obs_df,
      aes(x = dilutions, ymin = lo, ymax = hi),
      width = 0
    ) +
    geom_point(
      data = obs_df,
      aes(x = dilutions, y = survival),
      color = "black",
      size = 2
    ) +
    geom_point(
      data = data.frame(dilution = rm_endpoint, survival = 0.5),
      aes(x = dilution, y = survival),
      color = "orange2",
      fill = "orange2",
      shape = 23,
      size = 3.5
    ) +
    geom_point(
      data = data.frame(dilution = bayes_endpoint, survival = 0.5),
      aes(x = dilution, y = survival),
      color = "red3",
      shape = 8,
      size = 3.5
    ) +
    scale_x_log10(
      # breaks = raw_df$dilutions
      # ,labels = paste0("1:", raw_df$dilutions)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "Serum dilution",
      y = "Proportion surviving"
    ) +
    theme_minimal()
}

#### LOO-CV ####
make_loo_table <- function(loo_df) {
  gt::gt(loo_df) %>%
    gt::tab_header(
      title = "LOOCV comparison of SI and SII serocatalytic models"
    ) %>%
    gt::fmt_number(columns = c(ELPD_diff, SE_diff), decimals = 1) %>%
    gt::cols_label(
      Virus = "Enterovirus",
      Model = "Model Type",
      ELPD_diff = "ELPD diff.",
      SE_diff = "SE diff."
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "khaki1"),
      locations = gt::cells_body(
        rows = ELPD_diff == 0 & SE_diff == 0
      )
    )
}

#### Misc. ####

## SI Fig 1

plot_survival <- function(baseline_k1, baseline_phi) {
  
  # dilution seq
  d_seq <- seq(1, 350, by = 0.01)
  
  roma_colors <- color("roma")
  font_size <- 24
  
  k1_vals <- seq(0, 2, length.out=8)
  df_k1 <- expand.grid(d = d_seq, k1 = k1_vals)
  df_k1$survival_prob <- with(df_k1, survival_func(k1, baseline_phi, d))
  
  p1 <- ggplot(df_k1, aes(x = d, y = survival_prob, group = k1, color = k1)) +
    geom_line(lwd = 1.2) +
    scale_x_log10(breaks = c(1, 5, 10, 50, 100, 500, 1000)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.2)) +
    scale_color_gradientn(
      name = expression(k[1]),
      colours = roma_colors(length(k1_vals)),
      breaks = seq(0, 10, by = 2)
    ) +
    labs(x = "Dilution", y = "Survival") +
    theme_minimal() +
    theme(
      axis.line = element_line(colour = "black"),
      legend.key.height = unit(2, "cm"),
      legend.title.align = 0.5,
      legend.title = element_text(size = font_size, margin = margin(b = 30)),
      legend.text = element_text(size = font_size),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size)
    )
  
  phi_vals <- seq(1, 35, length.out = 8)
  df_phi <- expand.grid(d = d_seq, phi = phi_vals)
  df_phi$survival_prob <- with(df_phi, survival_func(baseline_k1, phi, d))
  
  p3 <- ggplot(df_phi, aes(x = d, y = survival_prob, group = phi, color = phi)) +
    geom_line(lwd = 1.2) +
    scale_x_log10(breaks = c(1, 5, 10, 50, 100, 500, 1000)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.2)) +
    scale_color_gradientn(
      name = expression(phi),
      colours = roma_colors(length(phi_vals)),
      breaks = seq(5, 35, by = 5)
    ) +
    labs(x = "Dilution", y = "Survival") +
    theme_minimal() +
    theme(
      axis.line = element_line(colour = "black"),
      legend.key.height = unit(2, "cm"),
      legend.title.align = 0.5,
      legend.title = element_text(size = font_size, margin = margin(b = 30)),
      legend.text = element_text(size = font_size),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size)
    )
  
  return(p3 / p1)
}

## SI Fig 2

plot_reed_muench_marginal_joint_bivariate <- function(reed_muench_fit) {
  font_size=20
  
  posterior_samples <- reed_muench_fit$draws(c("phi", "k1")) %>% as.data.frame()
  
  posterior_samples <- data.frame(phi = cbind(melt(posterior_samples[,1:4])$value),
                                  k1 = cbind(melt(posterior_samples[,5:8])$value))
  
  colnames(posterior_samples) <- c("phi", expression(k[1]))
  
  col_labels <- c(
    paste0(expression(phi), "\n"),
    expression(k[1])
  )
  
  ggpairs(
    posterior_samples,
    # columnLabels = col_labels,
    diag = list(continuous = wrap("barDiag", fill = "steelblue", bins=30)),
    upper = list(continuous = wrap("points", alpha = 0.5, col="steelblue")),
    lower = list(continuous = wrap("points", alpha = 0.5, col="steelblue")),
    labeller = label_parsed
  ) + theme_minimal() +
    theme(
      # panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.key.height = unit(2, "cm"),
      legend.title.align = 0.5,
      legend.title = element_text(size = font_size),
      legend.text = element_text(size = font_size),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      strip.text = element_text(size = font_size)
    )
  
}

plot_phi_vs_titer_combined <- function(phi_titer_summary_df) {
  
  df <- phi_titer_summary_df %>%
    dplyr::mutate(
      virus = factor(virus, levels = c("EV71", "CVA6", "EV68"))
    ) %>%
    dplyr::filter(titer < 2048) # for plotting
  
  df$virus <- factor(df$virus, levels = c("CVA6", "EV71", "EV68"))
  
  ggplot(df, aes(x = titer, color = virus)) +
    
    # vertical 95% credible intervals
    geom_linerange(
      aes(ymin = phi_lo, ymax = phi_hi),
      linewidth = 0.9
      ,position = position_dodge2(width=30)
    ) +
    
    # median points
    geom_point(
      aes(y = phi_med),
      size = 2.8
      ,position = position_dodge2(width=30)
    ) +
    
    scale_y_log10() +
    
    scale_x_continuous(
      breaks = c(0, 250, 500, 750, 1000),
      limits = c(0, 1100)
    ) +
    
    scale_virus_color(name = "Enterovirus") +
    
    labs(
      x = "Endpoint dilution titre",
      y = "Antibody concentration"
    ) +
    
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical"
    )
}

plot_phi_over_dilution <- function(df) {
  
  df <- df %>%
    dplyr::mutate(
      virus = factor(virus, levels = c("CVA6", "EV71", "EV68"))
    )
  
  ggplot(df, aes(x = dilutions, y = log10(phi_over_d), fill = virus, color = virus,
                 group = interaction(dilutions, virus))) +
    geom_boxplot(
      orientation = "x",
      position = position_dodge(width = 0.75),
      width = 0.55,
      outlier.shape = NA,
      linewidth = 0.6,
      color = "black"
    ) +
    scale_x_continuous(
      trans = "log2",
      breaks = sort(unique(df$dilutions))
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 6)
    ) +
    scale_virus_fill() +
    scale_virus_color() +
    labs(
      x = "Serum dilution",
      y = expression(log[10](phi[i] / d))
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)
    )
}

make_fig_2 <- function(reed_muench_plot, phi_vs_titer_plot, conc_vs_dilution_plot) {
  patchwork::wrap_plots(
    reed_muench_plot + 
      theme(
        axis.title = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)
      ),
    
    phi_vs_titer_plot + 
      theme(
        axis.title = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.position = c(0.8, 0.5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)
        ),
    
    conc_vs_dilution_plot + 
      theme(
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.text  = element_text(size = 14),
        legend.position = "none"
        ),
    
    nrow = 1,
    widths = c(1, 1, 1)
  ) +
    patchwork::plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
    theme(
      plot.tag = element_text(face = "bold")
    )
}

plot_phi_by_age_group <- function(phi_age_group_df) {
  
  age_levels <- c(
    "[0,5)",
    "[5,10)",
    "[10,20)",
    "[20,35)",
    "[35,55)",
    "[55,80)",
    "[80,95]"
  )
  
  df <- phi_age_group_df %>%
    dplyr::mutate(
      virus = factor(virus, levels = c("CVA6", "EV71", "EV68")),
      age_group = factor(age_group, levels = age_levels)
    )
  
  counts_df <- df %>%
    dplyr::group_by(virus, age_group) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  max_phi <- max(df$phi, na.rm = TRUE)
  label_x  <- min(max_phi * 1.10, 2700)
  
  ggplot(df, aes(x = phi, y = age_group, color = age_group)) +
    geom_boxplot(
      orientation = "y",
      outlier.shape = NA,
      width = 0.7,
      linewidth = 0.6,
      fill = NA
    ) +
    geom_text(
      data = counts_df,
      aes(x = label_x, y = age_group, label = paste0("n = ", n)),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = -0.6,   
      size = 2.5,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_viridis_d(option = "D", end = 0.9, guide = "none") +
    scale_x_continuous(
      # limits = c(0,2200)
    ) +
    facet_wrap(
      ~ virus,
      nrow = 1,
      labeller = as_labeller(virus_labels)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Antibody concentration",
      y = "Age (years)"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(
        fill = "grey85",
        color = "black",
        linewidth = 0.2
      ),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
      panel.spacing = grid::unit(0.8, "lines"),
      plot.margin = margin(5.5, 45, 5.5, 5.5)
    )
}



