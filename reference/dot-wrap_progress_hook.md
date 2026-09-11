# Signal progress without perturbing the worker RNG stream

\`progressor()\` (and condition signalling in general) draws random
condition identifiers inside the worker; restoring \`.Random.seed\`
around the call keeps the MCMC draws independent of the progress
reporting.

## Usage

``` r
.wrap_progress_hook(hook)
```
