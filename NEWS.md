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
