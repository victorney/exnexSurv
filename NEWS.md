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
