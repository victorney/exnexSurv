# Changelog

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
