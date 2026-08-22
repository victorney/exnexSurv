# Plot parameter traces from an exnex_surv fit

Generates one bayesplot traceplot per parameter.

## Usage

``` r
# S3 method for class 'exnex_surv'
plot(x, parameters = NULL, ask = interactive(), ...)
```

## Arguments

- x:

  A fitted \`exnex_surv\` object.

- parameters:

  Optional character vector of parameter names to plot.

- ask:

  Should R pause between plots? Defaults to \`interactive()\`.

- ...:

  Unused.

## Value

Invisibly returns \`x\`.
