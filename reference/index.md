# Package index

## Fitting

Fit an EXNEX survival model to right-censored data.

- [`exnex_surv()`](https://victorney.github.io/exnexSurv/reference/exnex_surv.md)
  : Fit Bayesian EXNEX Survival Models

## Inference and model selection

Posterior summaries, survival curves, model comparison and simulation.

- [`survival_curves()`](https://victorney.github.io/exnexSurv/reference/survival_curves.md)
  : Survival curves from an exnex_surv fit
- [`plot(`*`<survival_exnex>`*`)`](https://victorney.github.io/exnexSurv/reference/plot.survival_exnex.md)
  : Plot survival curves from an exnex_surv fit
- [`median_survival()`](https://victorney.github.io/exnexSurv/reference/median_survival.md)
  : Posterior median survival time
- [`rmst()`](https://victorney.github.io/exnexSurv/reference/rmst.md) :
  Restricted mean survival time (RMST) from an exnex_surv fit
- [`compute_waic()`](https://victorney.github.io/exnexSurv/reference/compute_waic.md)
  : Pointwise log-likelihood and WAIC from an exnex_surv fit
- [`compare_waic()`](https://victorney.github.io/exnexSurv/reference/compare_waic.md)
  : Compare multiple exnex_surv fits by WAIC
- [`probability_superiority()`](https://victorney.github.io/exnexSurv/reference/probability_superiority.md)
  : Posterior probability that one group beats another
- [`simulate_data()`](https://victorney.github.io/exnexSurv/reference/simulate_data.md)
  : Simulate basket-trial log-normal survival data

## Methods for fitted objects

Summarise, print, and plot fitted `exnex_surv` objects.

- [`summary(`*`<exnex_surv>`*`)`](https://victorney.github.io/exnexSurv/reference/summary.exnex_surv.md)
  : Summarize posterior draws from an exnex_surv fit
- [`print(`*`<exnex_surv>`*`)`](https://victorney.github.io/exnexSurv/reference/print.exnex_surv.md)
  : Print an exnex_surv fit
- [`plot(`*`<exnex_surv>`*`)`](https://victorney.github.io/exnexSurv/reference/plot.exnex_surv.md)
  : Plot parameter traces from an exnex_surv fit

## Internal functions

These functions are internal and not exported; use with caution.

- [`cpp_exnex_gibbs()`](https://victorney.github.io/exnexSurv/reference/cpp_exnex_gibbs.md)
  : Main Gibbs Sampler for EXNEX Survival Models
- [`exnex_surv_bridge()`](https://victorney.github.io/exnexSurv/reference/exnex_surv_bridge.md)
  : Bridge connecting hardhat processed data to the C++ Gibbs Sampler
- [`new_exnex_surv()`](https://victorney.github.io/exnexSurv/reference/new_exnex_surv.md)
  : Constructor for exnex_surv Objects
- [`.extract_event_vector()`](https://victorney.github.io/exnexSurv/reference/dot-extract_event_vector.md)
  : Extract event vector from outcomes
- [`.extract_time_vector()`](https://victorney.github.io/exnexSurv/reference/dot-extract_time_vector.md)
  : Extract time vector from outcomes
- [`.make_chain_seeds()`](https://victorney.github.io/exnexSurv/reference/dot-make_chain_seeds.md)
  : Derive one deterministic seed per chain (restores the global RNG
  state)
- [`.run_chains_parallel_exnex()`](https://victorney.github.io/exnexSurv/reference/dot-run_chains_parallel_exnex.md)
  : Run multiple chains in parallel via PSOCK workers
- [`.run_single_chain_exnex()`](https://victorney.github.io/exnexSurv/reference/dot-run_single_chain_exnex.md)
  : Run a single chain by calling the C++ kernel once
