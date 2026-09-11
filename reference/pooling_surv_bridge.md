# Bridge connecting hardhat processed data to the C++ Gibbs Sampler

Orchestrates data preparation, validation, and passes it to the C++ MCMC
engine.

## Usage

``` r
pooling_surv_bridge(
  processed,
  pooling,
  priors,
  iter,
  warmup,
  chains,
  parallel_chains,
  group_col,
  original_data,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- processed:

  A list produced by \`hardhat::mold()\`.

- pooling:

  Pooling mode; see
  [`pooling_surv`](https://victorney.github.io/exnexSurv/reference/pooling_surv.md)().

- priors:

  Optional named list of prior hyperparameters; see
  [`pooling_surv()`](https://victorney.github.io/exnexSurv/reference/pooling_surv.md).

- iter:

  Total number of MCMC iterations.

- warmup:

  Number of warmup iterations.

- chains:

  Number of chains to run.

- parallel_chains:

  Number of chains to run in parallel at the R level.

- group_col:

  Name of the original group column.

- original_data:

  The original data frame (before hardhat processing).

- seed:

  Optional seed used to deterministically derive one seed per chain.

## Value

An \`pooling_surv\` object.

## Details

The group variable (specified by group_col) is identified from the
original data. Additional variables are treated as covariates in the
linear predictor.
