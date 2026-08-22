# exnexSurv 1.2.0

* New vignette *The EXNEX Model, Priors, and Data Augmentation* presenting the model,
  the priors and their customization, the data augmentation mechanism, the Gibbs sampler,
  and a worked example with a shrinkage figure.
* Improved package-level documentation (overview, references) and expanded `exnex_surv()`
  documentation with runnable examples.
* Existing vignettes now use the exported `plot()` method instead of internal helpers.
* Expanded README and added this NEWS file.
* CRAN-readiness hardening: explicit NAMESPACE exports, `URL`/`BugReports` fields in
  DESCRIPTION, and a spelling word list.

# exnexSurv 1.1.0

* Basket-specific (vector-valued) priors for `p_mix`, `m_nex`, and `v_nex`, allowing one
  value per basket, with scalar fallback.
* `v_beta` is now configurable through `priors` (prior variance of the regression
  coefficients).
* Multiple chains with R-level parallel execution (`chains` and `parallel_chains`).
* Robustness and performance improvements to the C++ Gibbs sampler and its validation.
