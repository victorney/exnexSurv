# Deprecated alias for EXNEX fitting

\`exnex_surv()\` is deprecated in favour of \[pooling_surv()\]. It emits
a one-line deprecation message and forces \`pooling = "exnex"\`, so
existing calls continue to fit the EXNEX hierarchy unchanged.

## Usage

``` r
exnex_surv(x, ...)
```

## Arguments

- x:

  An object containing the predictors (subgroup assignment). Can be a
  data frame or a formula.

- ...:

  Additional arguments passed to methods.

## Value

An object of class \`pooling_surv\`, identical to \`pooling_surv(...,
pooling = "exnex")\`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- exnex_surv(survival::Surv(time, event) ~ group, data = d)
} # }
```
