library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gt)
library(matrixStats)

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
  dd <- draws_combined_melt %>% filter(variable %in% keep)
  
  var_labs <- c(
    "mu[1]" = "mu[S]",
    "mu[2]" = "mu[I[\"++\"]]",
    "mu[3]" = "mu[I[\"+\"]]",
    "sigma" = "sigma",
    "psi" = "psi",
    "lambda_1" = "lambda[1]",
    "beta" = "kappa"
  )
  
  ggplot(dd, aes(x = value, color = Virus, fill = Virus)) +
    geom_density(alpha = 0.25) +
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
      legend.position = "bottom",
      strip.text = element_text(size = 18, face = "bold")
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
      title = paste0("Parameter Summary: ", pretty_virus(virus))
    ) %>%
    fmt_markdown(columns = variable) %>%
    fmt_number(
      columns = c(mean, median, sd, q5, q95, rhat),
      decimals = 3
    ) %>%
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
    tab_style(
      style = cell_text(font = "Times New Roman", size = px(12)),
      locations = cells_body(everything())
    ) %>%
    tab_style(
      style = cell_text(font = "Times New Roman", size = px(13), weight = "bold"),
      locations = cells_title(groups = "title")
    ) %>%
    cols_width(
      variable ~ px(105),
      mean ~ px(68), median ~ px(68), sd ~ px(68),
      q5 ~ px(68), q95 ~ px(68), rhat ~ px(58),
      ess_bulk ~ px(72), ess_tail ~ px(72)
    ) %>%
    tab_options(
      table.font.names = "Times New Roman",
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
  
  ggplot(df_long, aes(x = Age)) +
    geom_ribbon(
      aes(ymin = lo, ymax = hi, fill = State),
      alpha = 0.18
    ) +
    geom_line(
      aes(y = mean, color = State),
      linewidth = 1.3
    ) +
    scale_serostate_color(state_names = state_names) +
    scale_serostate_fill(state_names = state_names) +
    labs(
      title = paste0("Serodynamics: ", pretty_virus(virus)),
      x = "Age (Years)",
      y = "Proportion",
      color = "Serostate",
      fill = "Serostate"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_weighted_avg_phi <- function(
    virus,
    plot_data,
    draws_df,
    summary_df,
    n_samples = 200,
    seed = 1,
    eps = 0.95
) {
  
  n_states <- if (any(grepl("^mu_x\\[3\\]", summary_df$variable))) 3 else 2
  
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
  
  state_labels <- serostate_levels[seq_len(n_states)]
  
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
      Ipp = rowMeans(I_mat)
    )
  }
  
  weighted_phi <- rowSums(
    probs_age * matrix(phi_state, nrow = nrow(probs_age), ncol = n_states, byrow = TRUE)
  )
  
  line_df <- data.frame(
    Age = tps,
    weighted_phi = weighted_phi
  )
  
  ggplot() +
    geom_point(
      data = pts,
      aes(x = age, y = phi, color = state, shape = conf),
      alpha = 0.85
    ) +
    geom_line(
      data = line_df,
      aes(x = Age, y = weighted_phi),
      linewidth = 1.4
    ) +
    scale_y_log10() +
    scale_serostate_color(state_names = state_labels) +
    scale_shape_manual(values = c(">95%" = 16, "<95%" = 0)) +
    labs(
      title = paste0("Serostate-weighted mean ", expression(phi), ": ", pretty_virus(virus)),
      x = "Age (Years)",
      y = expression(phi),
      color = "Serostate",
      shape = "Confidence"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.title.y = element_text(
        size = 16,
        angle = 0,
        vjust = 0.5,
        face = "bold"
      )
    )
}
