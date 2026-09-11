# Posterior rank probabilities for each group

For every posterior draw, the groups are ranked by a chosen survival
summary (median survival, survival at a fixed time, or RMST) and the
frequency of each rank is returned. Useful as a decision aid: the row of
the output for rank 1 gives the posterior probability that each group is
the best.

## Usage

``` r
rank_probabilities(
  fit,
  newdata = NULL,
  function_of = "median",
  times = NULL,
  tmax = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted \`pooling_surv\` object.

- newdata:

  Optional data frame with one row per group used to fix covariates. If
  \`NULL\`, covariates are set to zero.

- function_of:

  Character; one of \`"median"\`, \`"survival"\` (alias \`"S_t"\`), or
  \`"rmst"\`.

- times:

  Scalar time when \`function_of = "survival"\`.

- tmax:

  Positive horizon when \`function_of = "rmst"\`.

- ...:

  Unused.

## Value

A data frame with columns \`group\`, \`rank\`, and \`prob\`. Rank 1 is
the largest summary value (best group).

## Examples

``` r
set.seed(2719)
trial_data <- simulate_data(n = 90, K = 3, censoring_rate = 0.3)
fit <- pooling_surv(survival::Surv(time, event) ~ group,
                    data = trial_data,
                    iter = 300, warmup = 150, chains = 2)
rank_probabilities(fit)
#>   group rank        prob
#> 1     1    1 0.010000000
#> 2     2    1 0.053333333
#> 3     3    1 0.936666667
#> 4     1    2 0.233333333
#> 5     2    2 0.706666667
#> 6     3    2 0.060000000
#> 7     1    3 0.756666667
#> 8     2    3 0.240000000
#> 9     3    3 0.003333333
```
