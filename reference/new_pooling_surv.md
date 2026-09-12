# Constructor for pooling_surv Objects

Creates a formal S3 object containing posterior samples and metadata
from Bayesian pooling model fitting (EXNEX, complete pooling, or no
pooling).

## Usage

``` r
new_pooling_surv(
  draws,
  data,
  pooling,
  priors,
  resolved_priors = priors,
  iter,
  warmup,
  thin = 1L,
  chains,
  blueprint
)
```

## Arguments

- draws:

  A data frame containing posterior samples. Columns are: theta_1, ...,
  theta_K, beta_1, ..., beta_P (if P \> 0), sigma2.

- data:

  A list containing: - time: Observed follow-up times - event: Event
  indicators (0/1) - group: Group assignments - X: Covariate matrix (can
  be empty) - n: Total number of observations - n_groups: Number of
  groups - n_covariates: Number of covariates - cov_names: Covariate
  column names

- pooling:

  The fitted pooling variant: one of `"exnex"`, `"complete"`, or
  `"none"`.

- priors:

  A list of prior specifications used for fitting.

- resolved_priors:

  A named list of the prior hyperparameters actually used (defaults
  merged with any overrides supplied in `priors`).

- iter:

  Total MCMC iterations performed.

- warmup:

  Number of warmup iterations discarded.

- thin:

  Thinning interval applied to the post-warmup draws.

- chains:

  Number of chains run.

- blueprint:

  The hardhat blueprint for the original formula/data structure.

## Value

An \`pooling_surv\` object (S3 class).
