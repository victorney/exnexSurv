# exnexSurv

<!-- badges: start -->
[![R-CMD-check](https://github.com/victorney/exnexSurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/victorney/exnexSurv/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Fast Bayesian **EX**changeable–**N**on-**EX**changeable (EXNEX) survival analysis for
basket trials via Rcpp Gibbs sampling and data augmentation.

`exnexSurv` fits an EXNEX hierarchical model for right-censored log-normal survival data.
Each basket's log-location effect is either drawn from a shared exchangeable component
(borrowing strength across baskets) or from a basket-specific non-exchangeable prior,
selected by a latent indicator. Censored event times are imputed from their
truncated-Normal conditional distribution, which makes every full conditional conjugate and
the systematic Gibbs scan exact. The sampler is implemented in C++/RcppArmadillo and is
roughly 25&times; faster than a marginalized Stan implementation of the same model on the
companion paper's simulation study.

## Installation

```r
# development version from GitHub
# install.packages("pak")
pak::pak("victorney/exnexSurv")

# or
# remotes::install_github("victorney/exnexSurv")
```

## Quick example

```r
library(exnexSurv)
library(survival)

set.seed(1)
n <- 120
group <- factor(rep(1:3, each = 40))
x1 <- rnorm(n)
eta <- rep(c(1.1, 1.6, 2.0), each = 40) + 0.5 * x1
true_time <- exp(eta + rnorm(n, 0, 0.6))
cens_time <- runif(n, 2, 9)
d <- data.frame(
  time   = pmin(true_time, cens_time),
  event  = as.integer(true_time <= cens_time),
  group  = group,
  x1     = x1
)

fit <- exnex_surv(
  Surv(time, event) ~ group + x1,
  data = d,
  iter = 2000, warmup = 1000, chains = 2, parallel_chains = 2
)

summary(fit)
```

## Documentation

* **Vignettes** (built with `browseVignettes("exnexSurv")`):
  * [The EXNEX Model, Priors, and Data Augmentation](vignettes/model-and-methods.Rmd) — the model, the priors and how to customize them, the data augmentation mechanism, the Gibbs sampler, and a worked example.
  * [Getting Started](vignettes/getting-started.Rmd) — installing, fitting, and inspecting fits.
  * [Convergence and Diagnostics](vignettes/convergence-and-diagnostics.Rmd) — practical MCMC stability checks.
  * [Running Chains in Parallel](vignettes/parallel-chains.Rmd) — multiple chains and R-level parallelism.

## Methodology

The statistical methodology is described in:

> Ney, V. (2026). *Bayesian EXNEX survival models with data augmentation for basket
> trials.* Manuscript.

The model extends the robust exchangeability design of Neuenschwander et al. (2016),
*Pharmaceutical Statistics*, 15(2), 123–134.
