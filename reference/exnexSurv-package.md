# exnexSurv: Bayesian EXNEX Models for Survival Analysis in Basket Trials

Fits Bayesian Exchangeable–Non-Exchangeable (EXNEX) survival models for
right-censored log-normal data in basket trials using a fast,
data-augmented Gibbs sampler written in C++ (Rcpp/RcppArmadillo).

For patient \\i\\, \\\log T_i=\theta\_{g\[i\]}+X_i^{\mathsf T}\beta+
\varepsilon_i\\ with \\\varepsilon_i\sim\mathcal N(0,\sigma^2)\\. A
latent indicator \\Z_j\\ selects, for each basket \\j\\, whether its
effect is drawn from the exchangeable component \\\mathcal
N(\mu,\tau^2)\\ or from a basket-specific non-exchangeable prior
\\\mathcal N(m\_{0j},v\_{0j})\\. Censored event times are handled by
data augmentation (imputing truncated-Normal latent log-times), which
makes every full conditional conjugate and the systematic Gibbs scan
exact.

The main entry point is
[`exnex_surv()`](https://victorney.github.io/exnexSurv/reference/exnex_surv.md).
See the package vignette *The EXNEX Model, Priors, and Data
Augmentation* for a full presentation of the model, the priors, and the
sampler.

## References

Neuenschwander, B., Wandel, S., Roychoudhury, S., & Bailey, S. (2016).
Robust exchangeability designs for early phase clinical trials with
multiple strata. *Pharmaceutical Statistics*, 15(2), 123–134.

Tanner, M. A., & Wong, W. H. (1987). The calculation of posterior
distributions by data augmentation. *Journal of the American Statistical
Association*, 82(398), 528–540.

## See also

[`exnex_surv`](https://victorney.github.io/exnexSurv/reference/exnex_surv.md),
[`summary.exnex_surv`](https://victorney.github.io/exnexSurv/reference/summary.exnex_surv.md),
[`print.exnex_surv`](https://victorney.github.io/exnexSurv/reference/print.exnex_surv.md),
[`plot.exnex_surv`](https://victorney.github.io/exnexSurv/reference/plot.exnex_surv.md)

## Author

**Maintainer**: Victor Ney <victorney@ime.usp.br>

Authors:

- Victor Ney <victorney@ime.usp.br>
