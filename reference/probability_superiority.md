# Posterior probability that one group beats another

Estimates the posterior probability that a certain scalar summary of the
survival distribution for group \`a\` is greater than the same summary
for group \`b\`. Supported summaries are the median survival time, the
survival probability at a fixed time \`S(t_0)\`, and the restricted mean
survival time (RMST). The function works draw-by-draw on the posterior,
so the probability is computed on the joint posterior of the two groups.

## Usage

``` r
probability_superiority(
  fit,
  a,
  b,
  function_of = "median",
  newdata = NULL,
  times = NULL,
  tmax = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted \`pooling_surv\` object.

- a:

  An integer index (or group label) of the first group.

- b:

  An integer index (or group label) of the second group.

- function_of:

  Character; one of \`"median"\`, \`"survival"\` (alias \`"S_t"\`), or
  \`"rmst"\`. Selects which summary is compared.

- newdata:

  Optional data frame with one row per group, used to fix covariates. If
  \`NULL\`, covariates are set to zero.

- times:

  A scalar time at which to evaluate the survival function when
  \`function_of = "survival"\`.

- tmax:

  A positive horizon when \`function_of = "rmst"\`.

- ...:

  Unused.

## Value

A named list with elements \`prob\` (posterior probability that the
summary of group \`a\` exceeds that of group \`b\`), \`summary\` (which
summary was used), and \`groups\` (the two groups compared).
