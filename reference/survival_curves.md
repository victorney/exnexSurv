# Survival curves from an exnex_surv fit

Computes posterior survival curves \\S(t) = \Pr(T \> t)\\ from a fitted
\`exnex_surv\` model for a log-normal AFT specification: \$\$\log T_i =
\theta\_{g\[i\]} + X_i^\top\beta + \varepsilon,\qquad
\varepsilon\sim\mathcal N(0,\sigma^2).\$\$ For a fixed time \\t\\ and
linear predictor \\\eta\\, the survival probability of the log-normal
distribution is \$\$S(t) = 1 - \Phi\\\left(\frac{\log t -
\eta}{\sigma}\right),\$\$ where \\\sigma^2\\ is the residual variance.
Posterior draws of \\(\theta, \beta, \sigma^2)\\ are propagated through
this expression to obtain a full posterior distribution of \\S(t)\\ at
each time point.

## Usage

``` r
survival_curves(fit, newdata = NULL, times = NULL, level = 0.95, ...)
```

## Arguments

- fit:

  A fitted \`exnex_surv\` object.

- newdata:

  Optional data frame with columns matching the covariates of the model.
  If it contains a \`group\` column, that is used for the group index;
  otherwise the first group is used for all rows.

- times:

  Optional numeric vector of times at which to evaluate the curves. If
  \`NULL\`, a sensible grid is built from the observed follow-up times.

- level:

  Credible-interval level (default \`0.95\`).

- ...:

  Unused.

## Value

An object of class \`survival_exnex\` (also a \`data.frame\`) with long
format columns: \`time\`, \`median\`, \`lower\`, \`upper\`, and
\`group\`.

## Details

If \`newdata\` is supplied, each row is evaluated at its own group and
covariate values; otherwise the covariates are fixed at zero and the
first group is used (with a warning if more than one group exists).
