# Compare multiple pooling_surv fits by WAIC

Computes \[compute_waic()\] for each supplied fit and reports them in a
single \`data.frame\`, sorted by ascending WAIC.

## Usage

``` r
compare_waic(..., digit = 2)
```

## Arguments

- ...:

  Two or more fitted \`pooling_surv\` objects.

- digit:

  Number of decimal places for the reported statistics.

## Value

A \`data.frame\` with one row per model and columns \`model\`, \`waic\`,
\`se_elpd_waic\`, \`lpd\`, \`p_waic\`, \`elpd_waic\`.

## Details

Fits are labelled with the variable names used in the call when supplied
unnamed (e.g. \`fit_exnex\`, \`fit_complete\`), or by the names given to
\`...\`; when even those are unavailable, the fitted \`pooling\` variant
is used, with a final \`fit1\`/\`fit2\` fallback for duplicate or
anonymous labels.
