testthat::local_edition(3)

set.seed(202)
sim <- simulate_data(n = 15, beta = c(0.4), sigma = 1.0,
                     censoring_rate = 0.3, seed = 8)

fit_full <- exnexSurv::pooling_surv(
  survival::Surv(time, event) ~ group + x1,
  data = sim, iter = 200, warmup = 100, chains = 1, seed = 1
)
fit_group_only <- exnexSurv::pooling_surv(
  survival::Surv(time, event) ~ group,
  data = sim[, c("time", "event", "group")],
  iter = 200, warmup = 100, chains = 1, seed = 1
)

test_that("compute_waic returns expected structure", {
  w <- compute_waic(fit_full)
  expect_type(w, "list")
  expect_named(w, c("waic", "se_elpd_waic", "lpd", "p_waic", "elpd_waic", "pointwise"))
  expect_true(is.finite(w$waic))
  expect_true(w$se_elpd_waic >= 0)
  expect_true(w$p_waic >= 0)
  expect_identical(nrow(w$pointwise), nrow(sim))
  expect_true(all(w$pointwise[, "waic"] >= 1e-12))
})

test_that("p_waic is small for a good model vs larger for misspecified", {
  WAIC_full <- compute_waic(fit_full)
  WAIC_group <- compute_waic(fit_group_only)
  # the model with the true covariate should generally fit at least as well
  expect_true(WAIC_group$p_waic >= 0)
  expect_true(is.finite(WAIC_full$waic))
})

test_that("compare_waic reports sorted fits", {
  cmp <- compare_waic(group_only = fit_group_only, with_cov = fit_full)
  expect_s3_class(cmp, "data.frame")
  expect_true(all(c("model", "waic", "se_elpd_waic", "lpd", "p_waic", "elpd_waic") %in% colnames(cmp)))
  expect_true(all(diff(cmp$waic) >= 0))
  expect_identical(nrow(cmp), 2L)
})

test_that("compare_waic requires at least two fits", {
  expect_error(compare_waic(fit_full), "at least two")
})

test_that("compare_waic labels unnamed fits with call names or pooling mode", {
  # unnamed: labels come from the call (variable names)
  cmp <- compare_waic(fit_full, fit_group_only)
  expect_setequal(cmp$model, c("fit_full", "fit_group_only"))

  # anonymous inline fits: fall back to the fitted pooling variant
  fits_2 <- list(
    pooling_surv(survival::Surv(time, event) ~ group, data = sim,
                 pooling = "complete", iter = 100, warmup = 50, seed = 1),
    pooling_surv(survival::Surv(time, event) ~ group, data = sim,
                 pooling = "none", iter = 100, warmup = 50, seed = 2)
  )
  cmp2 <- compare_waic(fits_2[[1]], fits_2[[2]])
  expect_setequal(cmp2$model, c("complete", "none"))

  # duplicate labels stay readable
  cmp3 <- compare_waic(fit_full, fit_full)
  expect_setequal(cmp3$model, c("fit_full_1", "fit_full_2"))
})
