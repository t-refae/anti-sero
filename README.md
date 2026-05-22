# [Bayesian inference from virus neutralisation assay data]([url](https://doi.org/10.64898/2026.05.18.726027))

Bayesian pipeline from raw virus neutralisation assay data to age-dependent force of infection estimates for CVA6, EV-A71, and EV-D68.

A mechanistic antibody model estimates each individual's underlying antibody concentration (`φ`) from well-level survival across a dilution series. Posterior draws of `log φ` are then propagated (via Monte-Carlo integration in the Stan likelihood) into a serocatalytic model that infers age-varying FOI `λ(a) = λ₁exp(−κ(a−1))` and, for the three-state variant, a waning rate `ψ`.

## Models

- `Stan/antibody_mech.stan` — per-individual titers from binomial well counts, probit link, shared slope `k1`.
- `Stan/antibody_sero.stan` — **SII** (S → I++ → I+). Default for CVA6 and EV-A71.
- `Stan/antibody_sero_two_compartment.stan` — **SI** (S → I++). Default for EV-D68.

SI vs SII are compared per virus with approximate LOO-CV (`outputs/combined/combined_LOOCV_table.*`).

## Layout

```
_targets.R              # pipeline (targets + stantargets)
R/                      # 00_utils, 01_ppc, 02_data, 03_combine, 04_plots
Stan/                   # the three .stan files above
data/raw/               # {CA6,EV71,E68}_raw_well_observations.csv
outputs/                # per-virus + combined/ + reed_muench/ + misc/
```

Each raw CSV is one row per (individual × dilution) with `serumID`, `Age`, `Year`, `titer`, `rep1`, `rep2`, `outcome`, `dilutions`.

## Requirements

R ≥ 4.2, CmdStan (via `cmdstanr::install_cmdstan()`), a headless Chrome/Edge for `gt::gtsave()` PDF export, and the packages listed in `_targets.R`.

## Run

```r
targets::tar_make()
targets::tar_visnetwork()
```

OR

```r
source("run_and_time_pipeline.R")
```

MCMC settings and `n_phi_draws_thin` (the number of antibody posterior draws propagated into the sero stage) are at the top of `_targets.R`.
