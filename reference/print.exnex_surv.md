# Print an exnex_surv fit

Prints a compact summary and optionally traceplots of posterior draws.

## Usage

``` r
# S3 method for class 'exnex_surv'
print(x, show_trace = TRUE, parameters = NULL, max_parameters = Inf, ...)
```

## Arguments

- x:

  A fitted \`exnex_surv\` object.

- show_trace:

  Should traceplots be shown? Default is \`TRUE\`.

- parameters:

  Optional character vector of parameter names for traceplots.

- max_parameters:

  Maximum number of parameters to plot.

- ...:

  Unused.

## Value

Invisibly returns \`x\`.
