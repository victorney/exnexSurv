# Visualize shrinkage of group effects

Compares the hierarchical (EXNEX posterior) estimates of the
group-specific log-median survival \`theta_j\` (at zero covariates)
against the unpooled, basket-only Kaplan-Meier median estimate computed
from the original data. A short segment between the two points per group
shows the amount of shrinkage (partial pooling) applied by the model.

## Usage

``` r
shrinkage_plot(fit, ...)
```

## Arguments

- fit:

  A fitted \`pooling_surv\` object.

- ...:

  Unused.

## Value

Invisibly returns a data frame with columns \`group\`,
\`posterior_mean\`, \`posterior_lwr\`, \`posterior_upr\`, and
\`unpooled\` – the log-scale values used in the plot. Requires
\`ggplot2\`.

## Examples

``` r
set.seed(2719)
trial_data <- simulate_data(n = 90, K = 3, censoring_rate = 0.3)
fit <- pooling_surv(survival::Surv(time, event) ~ group,
                    data = trial_data,
                    iter = 300, warmup = 150, chains = 2)
if (FALSE) shrinkage_plot(fit) # \dontrun{}
```
