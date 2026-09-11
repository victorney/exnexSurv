# Restricted mean survival time (RMST) from an pooling_surv fit

Computes the restricted mean survival time up to a horizon \`tmax\`:
\$\$RMST(tmax) = \int_0^{tmax} S(t)\\dt.\$\$ The integral is evaluated
numerically (trapezoidal rule over a fine grid) for every posterior
draw, then summarised with posterior quantiles.

## Usage

``` r
rmst(fit, tmax = NULL, newdata = NULL, level = 0.95, grid_points = 400, ...)
```

## Arguments

- fit:

  A fitted \`pooling_surv\` object.

- tmax:

  A positive scalar horizon up to which the RMST is computed. Defaults
  to \`max(fit\$data\$time)\`.

- newdata:

  Optional data frame (one row gives one RMST).

- level:

  Credible-interval level (default \`0.95\`).

- grid_points:

  Resolution of the numerical quadrature (default \`400\`).

- ...:

  Unused.

## Value

A \`data.frame\` with columns \`group\`, \`rmst\`, \`lower\`, \`upper\`.
