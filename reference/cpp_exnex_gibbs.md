# Main Gibbs Sampler for EXNEX Survival Models

Implements Gibbs sampling for Bayesian EXNEX right-censored log-normal
survival models with optional covariates.

## Usage

``` r
cpp_exnex_gibbs(time, event, group, X, priors, iter, warmup, chains)
```

## Arguments

- time:

  Vector of observed follow-up times (n-vector, all \> 0)

- event:

  Vector of event indicators (n-vector, 0 or 1)

- group:

  Vector of group assignments (n-vector, integers 1 to K)

- X:

  Matrix of covariates (n x P). Can be empty (n x 0) if no covariates.

- priors:

  Optional named list of prior hyperparameters. Supported fields:
  `a_sigma`, `b_sigma`, `a_tau`, `b_tau` (inverse-Gamma shape and
  scale), `p_mix` (EXNEX mixture weight), `m_mu`, `v_mu` (exchangeable
  mean prior), `m_nex`, `v_nex` (nonexchangeable component), `v_beta`
  (variance of the regression-coefficient prior). `p_mix`, `m_nex`, and
  `v_nex` each accept either a scalar, replicated across baskets, or a
  numeric vector of length K with one value per basket, matching the
  basket-specific notation \\p\_{\mathrm{exch},j}\\, \\m\_{0j}\\,
  \\v\_{0j}\\ of the model. Absent fields keep the defaults; unknown
  fields are ignored.

- iter:

  Total number of MCMC iterations

- warmup:

  Number of iterations to discard

- chains:

  Number of independent chains to run

## Value

List containing posterior draws, priors, metadata, and diagnostics
