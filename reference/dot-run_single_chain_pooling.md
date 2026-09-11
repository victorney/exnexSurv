# Run a single chain by calling the C++ kernel once

\`rng_kind\` pins the RNG configuration (future workers may default to a
different normal-kind sampler, which would make the same seed produce
different draws than sequential execution).

## Usage

``` r
.run_single_chain_pooling(
  cpp_data,
  priors,
  pooling,
  iter,
  warmup,
  seed,
  verbose = TRUE,
  chain_label = "",
  progress_hook = NULL,
  rng_kind = NULL
)
```
