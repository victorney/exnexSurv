# Plot parameter traces from an pooling_surv fit

Generates one bayesplot traceplot per parameter.

## Usage

``` r
# S3 method for class 'pooling_surv'
plot(x, parameters = NULL, ask = interactive(), ...)
```

## Arguments

- x:

  A fitted \`pooling_surv\` object.

- parameters:

  Optional character vector of parameter names to plot.

- ask:

  Should R pause between plots? Defaults to \`interactive()\`.

- ...:

  Unused.

## Value

Invisibly returns \`x\`.
