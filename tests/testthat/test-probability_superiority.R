testthat::local_edition(3)

set.seed(303)
sim <- simulate_data(n_each = 20, beta = c(0.5), sigma = 1.1,
                     resistant = c(2, 8), resist_delta = -1.0,
                     censoring_rate = 0.3, seed = 12)
fit <- exnexSurv::exnex_surv(
  survival::Surv(time, event) ~ group + x1,
  data = sim, iter = 300, warmup = 150, chains = 1, seed = 2
)

test_that("probability_superiority on median is a valid probability", {
  p <- probability_superiority(fit, a = 1, b = 3, function_of = "median")
  expect_type(p, "list")
  expect_named(p, c("prob", "summary", "level", "groups"))
  expect_true(p$prob >= 0 & p$prob <= 1)
})

test_that("probability_superiority flips direction if arguments swapped", {
  p12 <- probability_superiority(fit, a = 1, b = 8, function_of = "median")
  p21 <- probability_superiority(fit, a = 8, b = 1, function_of = "median")
  expect_equal(p12$prob, 1 - p21$prob, tolerance = 1e-8)
})

test_that("probability_superiority supports survival at fixed time and rmst", {
  p_s <- probability_superiority(fit, 1, 2, function_of = "survival", times = 3)
  expect_true(p_s$prob >= 0 & p_s$prob <= 1)
  p_r <- probability_superiority(fit, 1, 2, function_of = "rmst", tmax = 10)
  expect_true(p_r$prob >= 0 & p_r$prob <= 1)
})

test_that("probability_superiority accepts group labels", {
  gn <- levels(sim$group)
  p_lab <- probability_superiority(fit, a = gn[1], b = gn[2], function_of = "median")
  p_num <- probability_superiority(fit, a = 1, b = 2, function_of = "median")
  expect_equal(p_lab$prob, p_num$prob)
})

test_that("probability_superiority validates args", {
  expect_error(probability_superiority(fit, 1, 1), "different groups")
  expect_error(probability_superiority(fit, 1, 2, function_of = "survival"),
               "`times` is required")
  expect_error(probability_superiority(fit, 1, 2, function_of = "rmst"),
               "`tmax` is required")
})
