# Package index

## Fitting

Fit a Bayesian pooling survival model: EXNEX (`pooling = "exnex"`),
complete pooling (`"complete"`), or no pooling (`"none"`).

- [`pooling_surv()`](https://victorney.github.io/exnexSurv/reference/pooling_surv.md)
  : Fit Bayesian Pooling Survival Models
- [`exnex_surv()`](https://victorney.github.io/exnexSurv/reference/exnex_surv.md)
  : Deprecated alias for EXNEX fitting

## Inference and model selection

Posterior summaries, survival curves, model comparison and simulation.

- [`survival_curves()`](https://victorney.github.io/exnexSurv/reference/survival_curves.md)
  : Survival curves from an pooling_surv fit
- [`plot(`*`<survival_exnex>`*`)`](https://victorney.github.io/exnexSurv/reference/plot.survival_exnex.md)
  : Plot survival curves from an pooling_surv fit
- [`median_survival()`](https://victorney.github.io/exnexSurv/reference/median_survival.md)
  : Posterior median survival time
- [`rmst()`](https://victorney.github.io/exnexSurv/reference/rmst.md) :
  Restricted mean survival time (RMST) from an pooling_surv fit
- [`compute_waic()`](https://victorney.github.io/exnexSurv/reference/compute_waic.md)
  : Pointwise log-likelihood and WAIC from an pooling_surv fit
- [`compare_waic()`](https://victorney.github.io/exnexSurv/reference/compare_waic.md)
  : Compare multiple pooling_surv fits by WAIC
- [`probability_superiority()`](https://victorney.github.io/exnexSurv/reference/probability_superiority.md)
  : Posterior probability that one group beats another
- [`rank_probabilities()`](https://victorney.github.io/exnexSurv/reference/rank_probabilities.md)
  : Posterior rank probabilities for each group
- [`convergence_diagnostics()`](https://victorney.github.io/exnexSurv/reference/convergence_diagnostics.md)
  : Convergence diagnostics (R-hat and effective sample size)
- [`simulate_data()`](https://victorney.github.io/exnexSurv/reference/simulate_data.md)
  : Simulate basket-trial log-normal survival data

## Model checking and visuals

Posterior predictive checks and shrinkage diagnostics for fitted
`pooling_surv` objects.

- [`posterior_predictive_check()`](https://victorney.github.io/exnexSurv/reference/posterior_predictive_check.md)
  : Posterior predictive survival curves
- [`shrinkage_plot()`](https://victorney.github.io/exnexSurv/reference/shrinkage_plot.md)
  : Visualize shrinkage of group effects

## Methods for fitted objects

Summarise, print, and plot fitted `pooling_surv` objects.

- [`summary(`*`<pooling_surv>`*`)`](https://victorney.github.io/exnexSurv/reference/summary.pooling_surv.md)
  : Summarize posterior draws from an pooling_surv fit
- [`print(`*`<pooling_surv>`*`)`](https://victorney.github.io/exnexSurv/reference/print.pooling_surv.md)
  : Print an pooling_surv fit
- [`plot(`*`<pooling_surv>`*`)`](https://victorney.github.io/exnexSurv/reference/plot.pooling_surv.md)
  : Plot parameter traces from an pooling_surv fit

## Internal functions

These functions are internal and not exported; use with caution.

- [`cpp_exnex_gibbs()`](https://victorney.github.io/exnexSurv/reference/cpp_exnex_gibbs.md)
  : Main Gibbs Sampler for EXNEX Survival Models
- [`pooling_surv_bridge()`](https://victorney.github.io/exnexSurv/reference/pooling_surv_bridge.md)
  : Bridge connecting hardhat processed data to the C++ Gibbs Sampler
- [`new_pooling_surv()`](https://victorney.github.io/exnexSurv/reference/new_pooling_surv.md)
  : Constructor for pooling_surv Objects
- [`.extract_event_vector()`](https://victorney.github.io/exnexSurv/reference/dot-extract_event_vector.md)
  : Extract event vector from outcomes
- [`.extract_time_vector()`](https://victorney.github.io/exnexSurv/reference/dot-extract_time_vector.md)
  : Extract time vector from outcomes
- [`.make_chain_seeds()`](https://victorney.github.io/exnexSurv/reference/dot-make_chain_seeds.md)
  : Derive one deterministic seed per chain (restores the global RNG
  state)
- [`.run_chains_parallel_pooling()`](https://victorney.github.io/exnexSurv/reference/dot-run_chains_parallel_pooling.md)
  : Run all chains concurrently on future workers
- [`.run_single_chain_pooling()`](https://victorney.github.io/exnexSurv/reference/dot-run_single_chain_pooling.md)
  : Run a single chain by calling the C++ kernel once
