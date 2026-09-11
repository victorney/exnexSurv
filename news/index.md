# Changelog

## exnexSurv 1.4.0

### New features

- [`pooling_surv()`](https://victorney.github.io/exnexSurv/reference/pooling_surv.md)
  is the new general entry point for fitting Bayesian pooling models for
  right-censored log-normal survival data. It runs the same
  data-augmented Gibbs sampler under three pooling modes, selected with
  the new `pooling` argument:
  - `pooling = "exnex"` — the EXNEX hierarchy (default, as in previous
    releases);
  - `pooling = "complete"` — a single shared basket effect with prior
    `N(m_mu, v_mu)` (complete pooling);
  - `pooling = "none"` — independent basket effects with priors
    `N(m_nex_j, v_nex_j)` (no pooling).

### Deprecation (soft)

- [`exnex_surv()`](https://victorney.github.io/exnexSurv/reference/exnex_surv.md)
  is deprecated in favour of
  [`pooling_surv()`](https://victorney.github.io/exnexSurv/reference/pooling_surv.md);
  it now emits a one-line deprecation message and forces
  `pooling = "exnex"`. All existing calls continue to work.
- Fitted objects are of class `"pooling_surv"`; all methods
  ([`summary()`](https://rdrr.io/r/base/summary.html),
  [`print()`](https://rdrr.io/r/base/print.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`survival_curves()`](https://victorney.github.io/exnexSurv/reference/survival_curves.md),
  [`median_survival()`](https://victorney.github.io/exnexSurv/reference/median_survival.md),
  [`rmst()`](https://victorney.github.io/exnexSurv/reference/rmst.md),
  [`probability_superiority()`](https://victorney.github.io/exnexSurv/reference/probability_superiority.md),
  [`compute_waic()`](https://victorney.github.io/exnexSurv/reference/compute_waic.md),
  [`compare_waic()`](https://victorney.github.io/exnexSurv/reference/compare_waic.md))
  accept them.

### Other

- Priors that are irrelevant for the chosen `pooling` mode (e.g. `p_mix`
  under `"complete"` or `"none"`) now trigger a warning naming the
  ignored fields.
- [`print()`](https://rdrr.io/r/base/print.html) output now labels the
  fitted pooling variant.
- New `verbose` argument in
  [`pooling_surv()`](https://victorney.github.io/exnexSurv/reference/pooling_surv.md)
  (default `TRUE`): a live progress bar (`progressr`) reports progress
  while the chains run; silenced with `verbose = FALSE`.
- [`compare_waic()`](https://victorney.github.io/exnexSurv/reference/compare_waic.md)
  now names rows with the argument names used in the call
  (e.g. `fit_exnex`) instead of the anonymous `fit1`/`fit2` labels; when
  no informative name is available it falls back to the fitted `pooling`
  variant.
- `parallel_chains` in
  [`pooling_surv()`](https://victorney.github.io/exnexSurv/reference/pooling_surv.md)
  is now a logical flag: `TRUE` runs all chains concurrently on
  background sessions through the `future` framework (progress from
  workers is relayed to your console, so it works in RStudio too);
  `FALSE` (default) runs chains sequentially.
- With `verbose = TRUE`, a live progress bar (`progressr`) now updates
  while the chains run, in both sequential and parallel execution; the
  per-chain text lines were removed.
- New
  [`convergence_diagnostics()`](https://victorney.github.io/exnexSurv/reference/convergence_diagnostics.md)
  computes split R-hat, bulk ESS and tail ESS per parameter;
  [`summary()`](https://rdrr.io/r/base/summary.html)/[`print()`](https://rdrr.io/r/base/print.html)
  now include these columns and a one-line convergence verdict
  automatically.
- New
  [`rank_probabilities()`](https://victorney.github.io/exnexSurv/reference/rank_probabilities.md)
  gives the posterior probability of each rank position per basket
  (median survival, `S(t_0)`, or RMST summaries).
- New
  [`posterior_predictive_check()`](https://victorney.github.io/exnexSurv/reference/posterior_predictive_check.md)
  overlays the posterior predictive mean survival curve with the
  observed Kaplan-Meier curve per basket.
- New
  [`shrinkage_plot()`](https://victorney.github.io/exnexSurv/reference/shrinkage_plot.md)
  contrasts EXNEX posterior estimates against unpooled Kaplan-Meier
  medians to visualise partial pooling.

## exnexSurv 1.3.1

### Breaking changes

- [`compute_waic()`](https://victorney.github.io/exnexSurv/reference/compute_waic.md)
  and
  [`compare_waic()`](https://victorney.github.io/exnexSurv/reference/compare_waic.md)
  now return `se_elpd_waic` instead of `se_waic`, matching the
  convention of the `loo` package.
- [`simulate_data()`](https://victorney.github.io/exnexSurv/reference/simulate_data.md)
  replaces argument `n_each` with `n` (scalar or per-basket vector) and
  adds argument `K` (number of baskets).

### Bug fixes

- [`simulate_data()`](https://victorney.github.io/exnexSurv/reference/simulate_data.md)
  now actually uses the `n` argument; previously `n_each` was accepted
  but ignored in favour of a hard-coded size sequence.
- [`survival_curves()`](https://victorney.github.io/exnexSurv/reference/survival_curves.md)
  now emits a warning when more than one group exists in the fitted
  model but `newdata` is not supplied.
- [`rmst()`](https://victorney.github.io/exnexSurv/reference/rmst.md)
  now correctly defaults `tmax` to `max(fit$data$time)` instead of
  requiring it.
- [`probability_superiority()`](https://victorney.github.io/exnexSurv/reference/probability_superiority.md)
  no longer returns a useless `level = NULL` element; its `newdata`
  handling is now consistent and better documented.

## exnexSurv 1.3.0

CRAN release: 2026-09-10

- New inference and model-selection tools (all in pure R, no changes to
  the validated C++ sampler):
  - [`survival_curves()`](https://victorney.github.io/exnexSurv/reference/survival_curves.md)
    computes posterior survival curves `S(t)` (with
    [`plot.survival_exnex()`](https://victorney.github.io/exnexSurv/reference/plot.survival_exnex.md)
    for plotting);
  - [`median_survival()`](https://victorney.github.io/exnexSurv/reference/median_survival.md)
    returns posterior median survival times;
  - [`rmst()`](https://victorney.github.io/exnexSurv/reference/rmst.md)
    returns the restricted mean survival time up to a horizon;
  - [`compute_waic()`](https://victorney.github.io/exnexSurv/reference/compute_waic.md)/[`compare_waic()`](https://victorney.github.io/exnexSurv/reference/compare_waic.md)
    provide WAIC model comparison;
  - [`probability_superiority()`](https://victorney.github.io/exnexSurv/reference/probability_superiority.md)
    estimates the posterior probability that one group beats another on
    median, survival probability, or RMST;
  - [`simulate_data()`](https://victorney.github.io/exnexSurv/reference/simulate_data.md)
    generates log-normal basket-trial datasets (with resistant baskets)
    for studies and examples.
- New vignette *Tools for Experimenters* covering these utilities.
