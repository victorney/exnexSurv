# Compare multiple exnex_surv fits by WAIC

Computes \[compute_waic()\] for each supplied fit and reports them in a
single \`data.frame\`, sorted by ascending WAIC.

## Usage

``` r
compare_waic(..., digit = 2)
```

## Arguments

- ...:

  Two or more fitted \`exnex_surv\` objects.

- digit:

  Number of decimal places for the reported statistics.

## Value

A \`data.frame\` with one row per model and columns \`model\`, \`waic\`,
\`se_elpd_waic\`, \`lpd\`, \`p_waic\`, \`elpd_waic\`.
