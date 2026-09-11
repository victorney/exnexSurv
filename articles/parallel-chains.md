# Running chains in parallel

## The `parallel_chains` flag

``` r

pooling_surv(..., chains = 4, parallel_chains = FALSE)  # sequential (default)
pooling_surv(..., chains = 4, parallel_chains = TRUE)   # concurrent
```

`parallel_chains` is a logical flag:

- `FALSE` (default): the `chains` chains run one after the other in the
  current R session.
- `TRUE`: all `chains` run concurrently on background R sessions via the
  [`future`](https://future.futureverse.org/) framework.

Changing `parallel_chains` affects only how the chains are scheduled,
not the model itself: the chains are independent and identical in both
modes, and the combination rule for the posterior draws is unchanged.

``` r

library(survival)
library(exnexSurv)

fit <- pooling_surv(
  survival::Surv(time, event) ~ group + covariate_1,
  data = trial_data,
  iter = 3000, warmup = 1500,
  chains = 4, parallel_chains = FALSE, seed = 42
)
```

## Progress reporting

With `verbose = TRUE` (default) a live progress bar (`progressr`)
updates while the chains run, in both sequential and parallel execution.

Set `verbose = FALSE` to run silently in both modes. Chains always run
with `iter` and `warmup` untouched by the scheduling mode, and `seed`
gives reproducibility in both modes via one seed per chain.

## Requirements

`parallel_chains = TRUE` uses the `future`, `future.apply`, and
`progressr` packages, which are installed automatically as dependencies
of `exnexSurv`.

The number of concurrent chains equals `chains`; there is no separate
worker-count option. If you want fewer concurrent runs, use
`parallel_chains = FALSE`.
