# Main Gibbs Sampler for EXNEX Survival Models

Implements Gibbs sampling for Bayesian EXNEX right-censored log-normal
survival models with optional covariates.

## Usage

``` r
cpp_exnex_gibbs(
  time,
  event,
  group,
  X,
  priors,
  pooling,
  verbose,
  iter,
  warmup,
  chains,
  chain_label,
  progress_hook = NULL
)
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

- pooling:

  Pooling mode: one of `"exnex"` (EXNEX hierarchy with mixture
  indicators \\Z_j\\, exchangeable center \\\mu\\ and between-basket
  variance \\\tau^2\\), `"complete"` (a single shared basket effect with
  prior \\\mathcal N(m\_{\mu}, v\_{\mu})\\; the EXNEX-specific priors
  are ignored) or `"none"` (independent basket effects with priors
  \\\mathcal N(m\_{0j}, v\_{0j})\\; the EXNEX-specific priors are
  ignored).

- verbose:

  If `TRUE`, report chain progress every 10% of the iterations on the
  standard output (or through `progress_hook`); silence with
  `verbose = FALSE`. Default `TRUE`.

- iter:

  Total number of MCMC iterations

- warmup:

  Number of iterations to discard

- chains:

  Number of independent chains to run

- chain_label:

  Optional label (e.g. the chain number) shown in the progress lines;
  empty string omits it.

- progress_hook:

  Optional R function called with one argument (the progress line) at
  every 10% milestone instead of printing to the standard output; used
  by the parallel runner to relay progress to the master session.

## Value

List containing posterior draws, priors, metadata, and diagnostics
