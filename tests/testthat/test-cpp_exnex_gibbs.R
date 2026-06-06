expect_error_message <- function(expr, pattern) {
  err <- tryCatch(expr, error = function(e) e)
  testthat::expect_s3_class(err, "error")
  testthat::expect_match(conditionMessage(err), pattern, fixed = TRUE)
}

test_that("cpp_exnex_gibbs returns expected draws and diagnostics", {
  res <- exnexSurv:::cpp_exnex_gibbs(
    time = c(5, 8, 12),
    event = c(1, 0, 1),
    group = c(1, 2, 1),
    X = matrix(nrow = 3, ncol = 0),
    priors = list(alpha = 1),
    iter = 4,
    warmup = 2,
    chains = 1
  )

  expect_type(res, "list")
  expect_identical(
    names(res),
    c("draws", "priors", "iter", "warmup", "chains", "diagnostics")
  )
  expect_true(is.matrix(res$draws))
  expect_identical(dim(res$draws), c(2L, 3L))
  expect_identical(colnames(res$draws), c("theta_1", "theta_2", "sigma"))
  expect_identical(
    unname(as.matrix(res$draws)),
    matrix(c(0, 0, 0, 0, 1, 1), nrow = 2)
  )
  expect_equal(res$diagnostics$n_obs, 3)
  expect_equal(res$diagnostics$n_groups, 2)
  expect_equal(res$diagnostics$n_covariates, 0)
  expect_equal(res$diagnostics$n_censored, 1)
  expect_equal(res$diagnostics$n_observed, 2)
  expect_identical(res$priors, list(alpha = 1))
})

test_that("cpp_exnex_gibbs validates input edge cases", {
  base_args <- list(
    time = c(5, 8, 12),
    event = c(1, 0, 1),
    group = c(1, 2, 1),
    X = matrix(c(60, 55, 62), ncol = 1),
    priors = list(),
    iter = 4,
    warmup = 2,
    chains = 1
  )

  expect_error_message(
    do.call(exnexSurv:::cpp_exnex_gibbs, modifyList(base_args, list(time = numeric(0)))),
    "time must have positive length."
  )

  expect_error_message(
    do.call(exnexSurv:::cpp_exnex_gibbs, modifyList(base_args, list(time = c(5, Inf, 12)))),
    "All survival times must be positive and finite."
  )

  expect_error_message(
    do.call(exnexSurv:::cpp_exnex_gibbs, modifyList(base_args, list(time = c(5, 0, 12)))),
    "All survival times must be positive and finite."
  )

  expect_error_message(
    do.call(exnexSurv:::cpp_exnex_gibbs, modifyList(base_args, list(iter = 0))),
    "iter must be positive."
  )

  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      modifyList(base_args, list(warmup = -1))
    ),
    "warmup must be non-negative."
  )

  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      modifyList(base_args, list(warmup = 4))
    ),
    "warmup must be strictly less than iter."
  )

  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      modifyList(base_args, list(event = c(1, 0)))
    ),
    "event must have the same length as time."
  )

  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      modifyList(base_args, list(event = c(1, 2, 0)))
    ),
    "event must contain only 0/1 values."
  )

  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      modifyList(base_args, list(group = c(1, 2)))
    ),
    "group must have the same length as time."
  )

  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      modifyList(base_args, list(X = matrix(c(60, 55), nrow = 2, ncol = 1)))
    ),
    "X must have the same number of rows as time."
  )

  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      modifyList(base_args, list(group = c(1, 2.5, 1)))
    ),
    "Group must contain integer values only."
  )

  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      modifyList(base_args, list(group = c(0, 2, 1)))
    ),
    "All group values must be >= 1."
  )
})
