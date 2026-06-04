test_that("default method rejects unsupported inputs", {
  expect_error(
    exnex_surv(1),
    "is not defined for a 'numeric'"
  )
})

test_that("formula interface returns a well-formed exnex_surv object without covariates", {
  df <- data.frame(
    group = c("A", "B", "A", "B"),
    time = c(1, 2, 3, 4),
    event = c(1, 0, 1, 0)
  )

  fit <- exnex_surv(
    survival::Surv(time, event) ~ group,
    data = df,
    iter = 10,
    warmup = 2,
    chains = 1
  )

  expect_s3_class(fit, "exnex_surv")
  expect_s3_class(fit, "hardhat_model")
  expect_named(
    fit,
    c("draws", "data", "priors", "iter", "warmup", "chains", "blueprint")
  )
  expect_s3_class(fit$blueprint, "hardhat_blueprint")
  expect_equal(fit$iter, 10)
  expect_equal(fit$warmup, 2)
  expect_equal(fit$chains, 1)
  expect_equal(fit$data$n_groups, 2)
  expect_equal(fit$data$n_covariates, 0)
  expect_equal(fit$data$cov_names, character(0))
  expect_equal(dim(fit$draws), c(8L, 3L))
  expect_named(fit$draws, c("theta_1", "theta_2", "sigma"))
})

test_that("formula interface keeps covariates and returns extra beta draws", {
  df <- data.frame(
    group = c("A", "B", "A", "B"),
    x1 = c(1, 2, 3, 4),
    time = c(1, 2, 3, 4),
    event = c(1, 0, 1, 0)
  )

  fit <- exnex_surv(
    survival::Surv(time, event) ~ group + x1,
    data = df,
    iter = 10,
    warmup = 2,
    chains = 1
  )

  expect_equal(fit$data$n_groups, 2)
  expect_equal(fit$data$n_covariates, 1)
  expect_equal(fit$data$cov_names, "x1")
  expect_equal(dim(fit$data$X), c(4L, 1L))
  expect_equal(dim(fit$draws), c(8L, 4L))
  expect_named(fit$draws, c("theta_1", "theta_2", "beta_1", "sigma"))
})

test_that("formula interface can use an explicit group_col that is not first in the formula", {
  df <- data.frame(
    x1 = c(10, 11, 12, 13),
    grp = c(1, 2, 1, 2),
    time = c(1, 2, 3, 4),
    event = c(1, 0, 1, 0)
  )

  fit <- exnex_surv(
    survival::Surv(time, event) ~ x1 + grp,
    data = df,
    group_col = "grp",
    iter = 10,
    warmup = 2,
    chains = 1
  )

  expect_equal(fit$data$group, c(1L, 2L, 1L, 2L))
  expect_equal(fit$data$n_groups, 2)
  expect_equal(fit$data$n_covariates, 1)
  expect_equal(fit$data$cov_names, "x1")
})

test_that("data.frame interface accepts Surv outcomes and returns the same structure", {
  x <- data.frame(
    group = c("A", "B", "A", "B"),
    x1 = c(1, 2, 3, 4)
  )
  y <- survival::Surv(c(1, 2, 3, 4), c(1, 0, 1, 0))

  fit <- exnex_surv(x, y, iter = 10, warmup = 2, chains = 1)

  expect_s3_class(fit, "exnex_surv")
  expect_equal(fit$data$n, 4)
  expect_equal(fit$data$n_groups, 2)
  expect_equal(fit$data$n_covariates, 1)
  expect_equal(fit$data$cov_names, "x1")
  expect_equal(dim(fit$draws), c(8L, 4L))
})

test_that("data.frame interface accepts matrix outcomes", {
  x <- data.frame(
    group = c("A", "B", "A", "B")
  )
  y <- cbind(time = c(1, 2, 3, 4), event = c(1, 0, 1, 0))

  fit <- exnex_surv(x, y, iter = 10, warmup = 2, chains = 1)

  expect_s3_class(fit, "exnex_surv")
  expect_equal(fit$data$n_groups, 2)
  expect_equal(fit$data$n_covariates, 0)
  expect_equal(dim(fit$draws), c(8L, 3L))
})

test_that("internal C++ wrapper returns the expected structure", {
  out <- exnexSurv:::cpp_exnex_gibbs(
    time = c(1, 2, 3, 4),
    event = c(1, 0, 1, 0),
    group = c(1, 2, 1, 2),
    X = matrix(numeric(0), nrow = 4, ncol = 0),
    priors = list(),
    iter = 10,
    warmup = 2,
    chains = 1
  )

  expect_named(
    out,
    c("draws", "priors", "iter", "warmup", "chains", "diagnostics")
  )
  expect_type(out$draws, "double")
  expect_equal(dim(out$draws), c(8L, 3L))
  expect_named(
    out$diagnostics,
    c("n_obs", "n_groups", "n_covariates", "n_censored", "n_observed")
  )
})

test_that("invalid inputs fail early with clear errors", {
  df <- data.frame(
    group = c("A", "B", "A", "B"),
    time = c(1, 2, 3, 4),
    event = c(1, 0, 1, 0)
  )

  expect_error(
    exnex_surv(
      survival::Surv(time, event) ~ group,
      data = df,
      iter = 10,
      warmup = 10,
      chains = 1
    ),
    "must be less than `iter`"
  )

  expect_error(
    exnex_surv(
      survival::Surv(time, event) ~ missing_group,
      data = df,
      iter = 10,
      warmup = 2,
      chains = 1
    ),
    "not found in data"
  )

  expect_error(
    exnex_surv(df, c(1, 0, 1), iter = 10, warmup = 2, chains = 1),
    "must be a Surv object, data frame, or matrix"
  )
})

test_that("formula interface errors when an explicit group_col is absent from data", {
  df <- data.frame(
    time = c(1, 2, 3, 4),
    event = c(1, 0, 1, 0)
  )

  expect_error(
    exnex_surv(survival::Surv(time, event) ~ x1, data = df, group_col = "missing_group", iter = 10, warmup = 2, chains = 1),
    "Column 'missing_group' not found in data"
  )
})
