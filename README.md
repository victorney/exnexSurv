# exnexSurv

<!-- badges: start -->
[![R-CMD-check](https://github.com/victorney/exnexSurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/victorney/exnexSurv/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Fast Bayesian logistic-pooling survival analysis for basket trials via Rcpp Gibbs
sampling and data augmentation. The `pooling` argument of `pooling_surv()` selects
between the **EX**changeable–**N**on-**EX**changeable (EXNEX) hierarchy
(`pooling = "exnex"`), complete pooling (`pooling = "complete"`), and no pooling
(`pooling = "none"`).

Under the EXNEX default, each basket's log-location effect is either drawn from a shared
exchangeable component (borrowing strength across baskets) or from a basket-specific
non-exchangeable prior, selected by a latent indicator. Censored event times are imputed from their
truncated-Normal conditional distribution, which makes every full conditional conjugate and
the systematic Gibbs scan exact. The sampler is implemented in C++/RcppArmadillo and is
roughly 25&times; faster than a marginalized Stan implementation of the same model on the 
paper's simulation study.

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

fit <- pooling_surv(
  Surv(time, event) ~ group + x1,
  data = d,
  iter = 2000, warmup = 1000, chains = 2, parallel_chains = 2
)

summary(fit)
```
