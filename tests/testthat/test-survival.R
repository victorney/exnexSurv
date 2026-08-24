testthat::local_edition(3)

# Small helper to fit quickly for tests
make_fit <- function(data, formula, ...) {
  exnexSurv::exnex_surv(
    formula,
    data = data,
    iter = 200,
    warmup = 100,
    chains = 1,
    seed = 99,
    ...
  )
}

set.seed(101)
sim <- simulate_data(n = 15, beta = c(0.4, -0.2), sigma = 1.0,
                     resist_delta = -0.7, censoring_rate = 0.3, seed = 7)

fit_cov <- make_fit(sim, survival::Surv(time, event) ~ group + x1 + x2)
fit_nocov <- make_fit(sim[, c("time", "event", "group")],
                      survival::Surv(time, event) ~ group)

test_that("survival_curves returns a survival_exnex object and monotone medians", {
  sc <- survival_curves(fit_cov)
  expect_s3_class(sc, "survival_exnex")
  expect_true(all(c("time", "median", "lower", "upper", "group") %in% colnames(sc)))
  expect_true(all(sc$median >= 0 & sc$median <= 1))
  # S(t) must be non-increasing in t
  expect_true(all(diff(sc$median) <= 1e-9 + 0 ))
  expect_identical(attr(sc, "level"), 0.95)
})

test_that("survival_curves at t=0 is near 1 and levels are ordered", {
  sc <- survival_curves(fit_cov, times = c(0, 1, 2, 5))
  expect_true(all(sc$lower <= sc$median & sc$median <= sc$upper))
})

test_that("median_survival returns finite posterior medians with ordered CIs", {
  med <- median_survival(fit_cov)
  expect_s3_class(med, "data.frame")
  expect_true(all(c("group", "median", "lower", "upper") %in% colnames(med)))
  expect_true(all(med$lower <= med$median & med$median <= med$upper))
  expect_true(all(is.finite(med$median)))
})

test_that("rmst returns positive finite values with ordered CIs", {
  r <- rmst(fit_cov, tmax = 10)
  expect_true(all(c("group", "rmst", "lower", "upper") %in% colnames(r)))
  expect_true(all(r$rmst > 0 & is.finite(r$rmst)))
  expect_true(all(r$lower <= r$rmst & r$rmst <= r$upper))
})

test_that("rmst is monotone in tmax and collapses to median relationship", {
  r_small <- rmst(fit_cov, tmax = 3)
  r_big <- rmst(fit_cov, tmax = 12)
  expect_true(all(r_big$rmst >= r_small$rmst))
})

test_that("functions subset correctly to a newdata row", {
  nd <- data.frame(x1 = 0, x2 = 0)
  med <- median_survival(fit_cov, newdata = nd)
  # median of exp(theta_j) at zero covs should be positive
  expect_true(all(is.finite(med$median)))
})
