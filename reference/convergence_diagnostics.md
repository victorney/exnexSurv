# Convergence diagnostics (R-hat and effective sample size)

Computes R-hat, bulk effective sample size (ESS), and tail ESS for every
parameter of a fitted \`pooling_surv\` object, using the \`posterior\`
package. All three follow the definitions in Vehtari, Gelman, Simpson,
Carpenter & Buerkner (2021, "Rank-normalization, folding, and
localization"): the R-hat is split and rank-normalized, bulk ESS is
computed on rank-normalized draws, and tail ESS quantifies sampling
adequacy in the 5 diagnostic is computed (so \`chains = 1\` still gets a
meaningful R-hat).

## Usage

``` r
convergence_diagnostics(fit, ...)
```

## Arguments

- fit:

  A fitted \`pooling_surv\` object (result of \`pooling_surv()\`).

- ...:

  Unused.

## Value

A data frame with one row per parameter and columns \`parameter\`,
\`rhat\`, \`ess_bulk\`, and \`ess_tail\`. All values are \`NA\` if the
\`posterior\` package is not installed.

## Details

Rules of thumb: R-hat close to 1.0 (say \< 1.01) and ESS in the hundreds
give reasonably stable summaries. Low ESS with a good R-hat means the
chain is slow to move but honest; both R-hat above ~1.01 and very low
ESS (e.g. below 100) are worth investigating before trusting the
posterior.

## Examples

``` r
set.seed(2719)
trial_data <- simulate_data(n = 90, K = 3, censoring_rate = 0.3)
fit <- pooling_surv(survival::Surv(time, event) ~ group,
                    data = trial_data,
                    iter = 300, warmup = 150, chains = 2)
convergence_diagnostics(fit)
#>   parameter      rhat ess_bulk ess_tail
#> 1   theta_1 1.0108117 195.3270 223.4138
#> 2   theta_2 1.0018604 136.8507 166.5135
#> 3   theta_3 0.9990224 185.2078 256.0480
#> 4    sigma2 1.0015918 150.8109 230.4038
```
