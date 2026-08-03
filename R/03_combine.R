# R/03_combine.R
library(dplyr)
library(tidyr)
library(reshape2)

combine_draws_three_viruses <- function(draws_ev71, draws_cva6, draws_ev68) {
  
  vars_3state <- c("mu[1]","mu[2]","mu[3]",
                   "sigma_phi[1]","sigma_phi[2]","sigma_phi[3]",
                   "psi","lambda_1","beta")
  vars_2state <- c("mu[1]","mu[2]",
                   "sigma_phi[1]","sigma_phi[2]",
                   "lambda_1","beta")
  
  ev71 <- draws_ev71 %>% dplyr::select(any_of(vars_3state)) %>% mutate(Virus="EV71")
  cva6 <- draws_cva6 %>% dplyr::select(any_of(vars_3state)) %>% mutate(Virus="CVA6")
  ev68 <- draws_ev68 %>% dplyr::select(any_of(vars_2state)) %>% mutate(Virus="EV68")
  
  bind_rows(ev71, cva6, ev68) %>%
    reshape2::melt(id.vars="Virus")
}

combine_summaries_three_viruses <- function(sum_ev71, sum_cva6, sum_ev68) {
  sum_ev71 <- sum_ev71 %>% mutate(Virus = "EV71")
  sum_cva6 <- sum_cva6 %>% mutate(Virus = "CVA6")
  sum_ev68 <- sum_ev68 %>% mutate(Virus = "EV68")
  bind_rows(sum_ev71, sum_cva6, sum_ev68)
}
