# Posterior predictive survival curves

Posterior predictive check of the fitted log-normal AFT model: the
posterior predictive mean survival curve per group (averaged over all
posterior draws) is overlaid with the observed Kaplan-Meier curve and a
band of replicated Kaplan-Meier curves. Systematic separation of the
observed KM from the posterior mean band indicates model misfit for that
group.

## Usage

``` r
posterior_predictive_check(
  fit,
  times = NULL,
  n_replicates = 200,
  newdata = NULL,
  plot = TRUE,
  ...
)
```

## Arguments

- fit:

  A fitted \`pooling_surv\` object.

- times:

  Optional vector of evaluation times; defaults to a grid of 100 points
  between 0 and the largest observed time.

- n_replicates:

  Number of posterior draws used for the replicated Kaplan-Meier band.
  Default 200.

- newdata:

  Optional data frame with one row per group giving the covariates for
  the predictive curves. If \`NULL\`, covariates are zero.

- plot:

  Should a \`ggplot2\` figure be drawn? Default \`TRUE\`.

- ...:

  Unused.

## Value

Invisibly returns a list with \`posterior\` (data frame: \`group\`,
\`time\`, \`mean\`, \`lwr\`, \`upr\`) and \`observed\` (data frame:
\`group\`, \`time\`, \`surv\`). Requires \`ggplot2\` when \`plot =
TRUE\`.

## Examples

``` r
set.seed(2719)
trial_data <- simulate_data(n = 90, K = 3, censoring_rate = 0.3)
fit <- pooling_surv(survival::Surv(time, event) ~ group,
                    data = trial_data,
                    iter = 300, warmup = 150, chains = 2)
if (FALSE) posterior_predictive_check(fit) # \dontrun{}
```
