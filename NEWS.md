# exnexSurv 1.4.0

## New features

* `pooling_surv()` is the new general entry point for fitting Bayesian
  pooling models for right-censored log-normal survival data. It runs the
  same data-augmented Gibbs sampler under three pooling modes, selected
  with the new `pooling` argument:
  * `pooling = "exnex"` — the EXNEX hierarchy (default, as in previous
    releases);
  * `pooling = "complete"` — a single shared basket effect with prior
    `N(m_mu, v_mu)` (complete pooling);
  * `pooling = "none"` — independent basket effects with priors
    `N(m_nex_j, v_nex_j)` (no pooling).

## Deprecation (soft)

* `exnex_surv()` is deprecated in favour of `pooling_surv()`; it now emits
  a one-line deprecation message and forces `pooling = "exnex"`. All
  existing calls continue to work.
* Fitted objects are of class `"pooling_surv"`; all methods
  (`summary()`, `print()`, `plot()`, `survival_curves()`,
  `median_survival()`, `rmst()`, `probability_superiority()`,
  `compute_waic()`, `compare_waic()`) accept them.

## Other

* Priors that are irrelevant for the chosen `pooling` mode (e.g. `p_mix`
  under `"complete"` or `"none"`) now trigger a warning naming the ignored
  fields.
* `print()` output now labels the fitted pooling variant.

# exnexSurv 1.3.1

## Breaking changes

* `compute_waic()` and `compare_waic()` now return `se_elpd_waic` instead
  of `se_waic`, matching the convention of the `loo` package.
* `simulate_data()` replaces argument `n_each` with `n` (scalar or
  per-basket vector) and adds argument `K` (number of baskets).

## Bug fixes

* `simulate_data()` now actually uses the `n` argument; previously `n_each`
  was accepted but ignored in favour of a hard-coded size sequence.
* `survival_curves()` now emits a warning when more than one group exists
  in the fitted model but `newdata` is not supplied.
* `rmst()` now correctly defaults `tmax` to `max(fit$data$time)` instead
  of requiring it.
* `probability_superiority()` no longer returns a useless `level = NULL`
  element; its `newdata` handling is now consistent and better documented.

# exnexSurv 1.3.0

* New inference and model-selection tools (all in pure R, no changes to the
  validated C++ sampler):
  * `survival_curves()` computes posterior survival curves `S(t)` (with
    `plot.survival_exnex()` for plotting);
  * `median_survival()` returns posterior median survival times;
  * `rmst()` returns the restricted mean survival time up to a horizon;
  * `compute_waic()`/`compare_waic()` provide WAIC model comparison;
  * `probability_superiority()` estimates the posterior probability that one
    group beats another on median, survival probability, or RMST;
  * `simulate_data()` generates log-normal basket-trial datasets (with
    resistant baskets) for studies and examples.
* New vignette *Tools for Experimenters* covering these utilities.
