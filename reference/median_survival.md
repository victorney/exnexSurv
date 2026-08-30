# Posterior median survival time

For the log-normal AFT model the median survival time for a linear
predictor \\\eta\\ is \\t\_{med} = \exp(\eta)\\, since \\S(t)=0.5\\ when
\\\log t = \eta\\. Posterior draws of \\\theta\\ and \\\beta\\ therefore
induce a posterior distribution of \\t\_{med}\\ whose quantiles are
reported.

## Usage

``` r
median_survival(fit, newdata = NULL, level = 0.95, ...)
```

## Arguments

- fit:

  A fitted \`exnex_surv\` object.

- newdata:

  Optional data frame (one row gives one median).

- level:

  Credible-interval level (default \`0.95\`).

- ...:

  Unused.

## Value

A \`data.frame\` with columns \`group\`, \`median\`, \`lower\`,
\`upper\`.
