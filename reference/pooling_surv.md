# Fit Bayesian Pooling Survival Models

Fits Bayesian models for right-censored log-normal survival data in
basket trials using a data-augmented Gibbs sampler. The \`pooling\`
argument selects the structural assumption on the basket effects:

## Usage

``` r
pooling_surv(x, ...)

# Default S3 method
pooling_surv(x, ...)

# S3 method for class 'formula'
pooling_surv(
  formula,
  data,
  pooling = "exnex",
  priors = list(),
  iter = 2000,
  warmup = 1000,
  chains = 1,
  parallel_chains = FALSE,
  group_col = NULL,
  seed = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'data.frame'
pooling_surv(
  x,
  y,
  pooling = "exnex",
  priors = list(),
  iter = 2000,
  warmup = 1000,
  chains = 1,
  parallel_chains = FALSE,
  group_col = NULL,
  seed = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object containing the predictors (subgroup assignment). Can be a
  data frame or a formula.

- ...:

  Additional arguments passed to methods.

- formula:

  A model formula with structure: \`Surv(time, event) ~ group +
  covariates\`. The first RHS variable is the group/basket assignment.
  Additional variables are covariates.

- data:

  A data frame containing the variables in the formula.

- pooling:

  Pooling mode for the basket effects: one of `"exnex"` (default),
  `"complete"`, or `"none"` (see Description).

- priors:

  Optional named list of prior hyperparameters. Supported fields:
  `a_sigma`, `b_sigma` (inverse-Gamma shape and scale for the residual
  variance \\\sigma^2\\); `a_tau`, `b_tau` (inverse-Gamma for the
  between-basket variance \\\tau^2\\. Only used with
  `pooling = "exnex"`); `p_mix` (EXNEX mixture weight, strictly between
  0 and 1. Only used with `pooling = "exnex"`); `m_mu`, `v_mu` (prior
  mean and variance of the exchangeable center \\\mu\\; with
  `pooling = "complete"` they act as the prior of the shared basket
  effect); `m_nex`, `v_nex` (prior mean and variance of the
  nonexchangeable component; with `pooling = "none"` they are the priors
  of the independent basket effects); `v_beta` (variance of the
  regression-coefficient prior). `p_mix`, `m_nex`, and `v_nex` each
  accept either a scalar, replicated across baskets, or a numeric vector
  of length K with one value per basket, matching the basket-specific
  notation of the model. All other fields are scalars. Absent fields
  keep the defaults (inverse-Gamma(2,2) for the variances,
  `p_mix = 0.5`, `m_mu = 0`, `v_mu = 1e4`, `m_nex = 0`, `v_nex = 1e4`,
  `v_beta = 1e4`); unknown fields are ignored. Priors that are
  irrelevant for the selected `pooling` mode are ignored with a warning.

- iter:

  Total number of MCMC iterations. Default is 2000.

- warmup:

  Number of warmup iterations to discard. Default is 1000. Posterior
  samples will have (iter - warmup) rows.

- chains:

  Number of independent MCMC chains. Default is 1.

- parallel_chains:

  Logical. If \`TRUE\`, all \`chains\` run concurrently on background R
  sessions (via the \`future\` framework); if \`FALSE\` (default),
  chains run sequentially, one after the other.

- group_col:

  Name of the column in \`x\` that represents the basket/group
  assignment. If NULL (default), assumes the first column in x is the
  group.

- seed:

  Random seed for reproducibility (optional).

- verbose:

  If TRUE (default), show a live progress bar (\`progressr\`) that
  advances as the chains run, in both sequential and parallel execution.
  Silence it with \`verbose = FALSE\`.

- y:

  A Surv object or matrix containing outcome (time and event status).

## Value

An object of class \`pooling_surv\` containing posterior samples and
metadata. Components include `pooling` (the fitted variant), `draws` (a
data frame with columns `theta_1`, ..., `theta_K`, `beta_1`, ...,
`beta_P`, and `sigma2`; under `pooling = "complete"` the `theta_*`
columns hold identical draws from the shared effect), `data` (the
processed data: `time`, `event`, `group`, `X`, `n`, `n_groups`,
`n_covariates`, `cov_names`, `chain_seeds`), `priors` (the supplied
priors), `resolved_priors` (defaults merged with overrides), `iter`,
`warmup`, `chains`, and `blueprint`. Use
[`summary()`](https://rdrr.io/r/base/summary.html),
[`print()`](https://rdrr.io/r/base/print.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) to inspect it.

## Details

For patient \\i\\, the model is \$\$\log T_i = \theta\_{g\[i\]} +
X_i^{\mathsf T}\beta + \varepsilon_i,\quad \varepsilon_i\sim\mathcal
N(0,\sigma^2),\$\$ and the argument `pooling` selects the prior
structure of the basket effects \\\theta_j\\:

- `"exnex"` (default):

  latent EXNEX hierarchy, \$\$\theta_j\mid Z_j\sim Z_j\\\mathcal
  N(\mu,\tau^2) +(1-Z_j)\\\mathcal N(m\_{0j},v\_{0j}),\quad
  Z_j\sim\mathrm{Bern}(p\_{\mathrm{exch},j}),\$\$ so that each basket
  either borrows strength from the exchangeable component \\\mathcal
  N(\mu,\tau^2)\\ or follows its own non-exchangeable prior \\\mathcal
  N(m\_{0j},v\_{0j})\\.

- `"complete"`:

  complete pooling: a single shared basket effect \\\theta_j=\theta\\
  for all baskets with prior \\\theta\sim\mathcal
  N(m\_{\mu},v\_{\mu})\\.

- `"none"`:

  no pooling: independent basket effects, \\\theta_j\sim\mathcal
  N(m\_{0j},v\_{0j})\\ for \\j=1,\dots,K\\.

The first predictor variable in the model formula is interpreted as the
basket/group assignment; additional variables are covariates in the
linear predictor on the log-survival scale.

Censored event times are handled by data augmentation: a censored
log-time is imputed from its truncated-Normal conditional distribution
(inverse-CDF in log space) before the remaining parameters are updated
with their conjugate full conditionals. The sampler keeps only the
posterior draws of \\\theta_j\\, \\\beta\\, and \\\sigma^2\\.

## Examples

``` r
# Small simulated dataset: three baskets, one covariate, ~25% censoring
set.seed(1)
n <- 90
group <- factor(rep(1:3, each = 30))
x1 <- rnorm(n)
eta <- rep(c(1.1, 1.6, 2.0), each = 30) + 0.5 * x1
true_time <- exp(eta + rnorm(n, 0, 0.6))
cens_time <- runif(n, 2, 9)
d <- data.frame(
  time = pmin(true_time, cens_time),
  event = as.integer(true_time <= cens_time),
  group = group,
  x1 = x1
)

# Fit the EXNEX hierarchy with group effects and one covariate
fit <- pooling_surv(
  survival::Surv(time, event) ~ group + x1,
  data = d,
  priors = list(p_mix = 0.7),
  iter = 300, warmup = 150, chains = 1
)
print(fit, show_trace = FALSE)
#> <pooling_surv model>
#> Pooling: exnex 
#> Draws: 150 total post-warmup samples
#>        150 post-warmup samples per chain
#> Groups: 3 | Covariates: 1 
#> MCMC: iter = 300 , warmup = 150 , chains = 1 
#> 
#>  parameter      mean         sd       q05       q50       q95     rhat ess_bulk
#>    theta_1 1.2335895 0.13974939 1.0095657 1.2187444 1.4675805 1.016317 78.54732
#>    theta_2 1.4724191 0.13478558 1.2427548 1.4719397 1.6794470 1.002595 70.19314
#>    theta_3 2.1979588 0.15277914 1.9438532 2.1936782 2.4377354 1.003304 56.67331
#>     beta_1 0.5189845 0.08791554 0.3809179 0.5196996 0.6528267 1.016858 92.05982
#>     sigma2 0.4029807 0.08911116 0.2882274 0.3903261 0.5883140 1.014836 30.91370
#>   ess_tail
#>  104.18341
#>   89.79518
#>  131.31411
#>  132.80604
#>   31.09511
#> 
#> Convergence: max R-hat = 1.02 | min ESS = 30.9 (consider longer chains or more warmup)
summary(fit)
#>   parameter      mean         sd       q05       q50       q95     rhat
#> 1   theta_1 1.2335895 0.13974939 1.0095657 1.2187444 1.4675805 1.016317
#> 2   theta_2 1.4724191 0.13478558 1.2427548 1.4719397 1.6794470 1.002595
#> 3   theta_3 2.1979588 0.15277914 1.9438532 2.1936782 2.4377354 1.003304
#> 4    beta_1 0.5189845 0.08791554 0.3809179 0.5196996 0.6528267 1.016858
#> 5    sigma2 0.4029807 0.08911116 0.2882274 0.3903261 0.5883140 1.014836
#>   ess_bulk  ess_tail
#> 1 78.54732 104.18341
#> 2 70.19314  89.79518
#> 3 56.67331 131.31411
#> 4 92.05982 132.80604
#> 5 30.91370  31.09511

# Complete and no pooling use the same data augmentation machinery
fit_complete <- pooling_surv(
  survival::Surv(time, event) ~ group + x1,
  data = d, pooling = "complete",
  iter = 300, warmup = 150, chains = 1
)
fit_none <- pooling_surv(
  survival::Surv(time, event) ~ group + x1,
  data = d, pooling = "none",
  iter = 300, warmup = 150, chains = 1
)
```
