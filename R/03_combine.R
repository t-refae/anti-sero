# R/03_combine.R
library(dplyr)
library(tidyr)
library(reshape2)

combine_draws_three_viruses <- function(draws_ev71, draws_cva6, draws_ev68) {
  
  vars <- c("mu_x[1]", "mu_x[2]", "mu_x[3]",
            "sigma_phi[1]", "sigma_phi[2]", "sigma_phi[3]",
            "psi", "lambda_1", "kappa")
  
  pick <- function(d, virus) {
    d %>%
      dplyr::select(dplyr::any_of(vars)) %>%
      dplyr::mutate(Virus = virus)
  }
  
  dplyr::bind_rows(
    pick(draws_ev71, "EV71"),
    pick(draws_cva6, "CVA6"),
    pick(draws_ev68, "EV68")
  ) %>%
    reshape2::melt(id.vars = "Virus", na.rm = TRUE)
}

combine_summaries_three_viruses <- function(sum_ev71, sum_cva6, sum_ev68) {
  sum_ev71 <- sum_ev71 %>% mutate(Virus = "EV71")
  sum_cva6 <- sum_cva6 %>% mutate(Virus = "CVA6")
  sum_ev68 <- sum_ev68 %>% mutate(Virus = "EV68")
  bind_rows(sum_ev71, sum_cva6, sum_ev68)
}
