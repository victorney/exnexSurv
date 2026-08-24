testthat::local_edition(3)

test_that("simulate_data returns expected structure and attributes", {
  set.seed(1)
  d <- simulate_data(n = 12, beta = c(0.5, -0.2), sigma = 1.0,
                     outlier_baskets = c(2, 5), resist_delta = -0.8,
                     censoring_rate = 0.3, seed = 5)
  expect_s3_class(d, "data.frame")
  expect_true(all(c("time", "event", "group", "x1", "x2") %in% colnames(d)))
  expect_true(all(d$time > 0))
  expect_true(all(d$event %in% c(0, 1)))
  expect_identical(attr(d, "true_beta"), c(0.5, -0.2))
  expect_equal(attr(d, "true_sigma"), 1.0)
  expect_length(attr(d, "true_theta"), 9)
})

test_that("simulate_data applies resist_delta to selected baskets", {
  set.seed(2)
  d <- simulate_data(n = 40, beta = c(0.5), sigma = 0.01,
                     outlier_baskets = 3, resist_delta = -1.0, seed = 6)
  tt <- attr(d, "true_theta")
  # basket-3 theta should be ~1 unit below the healthy population
  healthy <- mean(tt[-3])
  expect_equal(tt[3], healthy - 1.0, tolerance = 0.4)
})

test_that("simulate_data reproduces with the same seed", {
  d1 <- simulate_data(seed = 123)
  d2 <- simulate_data(seed = 123)
  expect_identical(d1, d2)
})

test_that("simulate_data validates inputs", {
  expect_error(simulate_data(n = 0))
  expect_error(simulate_data(sigma = -1), "sigma")
  expect_error(simulate_data(outlier_baskets = 100), "outlier_baskets")
})

test_that("simulate_data accepts n as a scalar with custom K", {
  d <- simulate_data(n = 25, K = 5, beta = c(0.3), seed = 9)
  expect_length(attr(d, "true_theta"), 5)
  expect_equal(nrow(d), 5 * 25)
})

test_that("simulate_data accepts n as a per-basket vector", {
  d <- simulate_data(n = c(10, 15, 20, 10, 12), K = 5, beta = c(0.3), seed = 10)
  expect_length(attr(d, "true_theta"), 5)
  expect_equal(nrow(d), sum(c(10, 15, 20, 10, 12)))
})

test_that("simulate_data validates n/K length mismatch", {
  expect_error(simulate_data(n = c(10, 20), K = 5),
               "must be a scalar")
})
