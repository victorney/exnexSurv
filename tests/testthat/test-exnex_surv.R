testthat::local_edition(3)

expect_error_message <- function(expr, pattern) {
  err <- tryCatch(expr, error = function(e) e)
  testthat::expect_s3_class(err, "error")
  testthat::expect_match(conditionMessage(err), pattern, fixed = TRUE)
}

trial_data <- data.frame(
  time = c(5, 8, 12, 9),
  event = c(1, 0, 1, 1),
  group = factor(c("A", "B", "A", "B")),
  age = c(60, 55, 62, 58)
)

expected_draws <- data.frame(
  theta_1 = c(0, 0),
  theta_2 = c(0, 0),
  beta_1 = c(0, 0),
  sigma = c(1, 1)
)

test_that("exnex_surv works with the formula interface and hardhat blueprints", {
  fit <- exnex_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    priors = list(alpha = 1),
    iter = 4,
    warmup = 2,
    chains = 1
  )

  expect_identical(
    class(fit),
    c("exnex_surv", "hardhat_model", "hardhat_scalar")
  )
  expect_identical(
    class(fit$blueprint),
    c("default_formula_blueprint", "formula_blueprint", "hardhat_blueprint")
  )
  expect_identical(fit$draws, expected_draws)
  expect_identical(fit$data$time, trial_data$time)
  expect_identical(fit$data$event, trial_data$event)
  expect_identical(fit$data$group, c(1, 2, 1, 2))
  expect_identical(fit$data$X[, 1], trial_data$age)
  expect_identical(fit$data$cov_names, "age")
  expect_equal(fit$data$n, 4)
  expect_equal(fit$data$n_groups, 2)
  expect_equal(fit$data$n_covariates, 1)
})

test_that("exnex_surv works with the data.frame interface", {
  fit <- exnex_surv(
    x = trial_data[c("group", "age")],
    y = survival::Surv(trial_data$time, trial_data$event),
    priors = list(alpha = 1),
    iter = 4,
    warmup = 2,
    chains = 1
  )

  expect_identical(
    class(fit),
    c("exnex_surv", "hardhat_model", "hardhat_scalar")
  )
  expect_identical(
    class(fit$blueprint),
    c("default_xy_blueprint", "xy_blueprint", "hardhat_blueprint")
  )
  expect_identical(fit$draws, expected_draws)
  expect_identical(fit$data$group, c(1, 2, 1, 2))
  expect_identical(fit$data$X[, 1], trial_data$age)
  expect_identical(fit$data$cov_names, "age")
})

test_that("exnex_surv keeps covariates that only share the group prefix", {
  prefixed_data <- data.frame(
    time = c(5, 8, 12, 9),
    event = c(1, 0, 1, 1),
    group = factor(c("A", "B", "A", "B")),
    group_size = c(10, 11, 12, 13),
    age = c(60, 55, 62, 58)
  )

  fit <- exnex_surv(
    survival::Surv(time, event) ~ group + group_size + age,
    data = prefixed_data,
    priors = list(alpha = 1),
    iter = 4,
    warmup = 2,
    chains = 1
  )

  expect_identical(colnames(fit$data$X), c("group_size", "age"))
  expect_identical(fit$data$cov_names, c("group_size", "age"))
  expect_identical(unname(fit$data$X[, "group_size"]), prefixed_data$group_size)
  expect_identical(unname(fit$data$X[, "age"]), prefixed_data$age)
  expect_equal(fit$data$n_covariates, 2)
})

test_that("exnex_surv validates public API edge cases", {
  expect_error_message(
    exnex_surv(
      survival::Surv(time, event) ~ group,
      data = trial_data,
      priors = list(),
      iter = 4,
      warmup = 4
    ),
    "`warmup` (4) must be less than `iter` (4)."
  )

  expect_error_message(
    exnex_surv(
      survival::Surv(time, event) ~ group,
      data = trial_data,
      priors = list(),
      iter = 4,
      warmup = 2,
      chains = 2
    ),
    "Multi-chain support is not implemented yet. Use chains = 1."
  )

  expect_error_message(
    exnex_surv(
      survival::Surv(time, event) ~ group,
      data = trial_data,
      priors = list(),
      iter = 4,
      warmup = 2,
      group_col = "missing"
    ),
    "Column 'missing' not found in data."
  )

  expect_error_message(
    exnex_surv(
      x = trial_data[c("group", "age")],
      y = trial_data$time,
      priors = list(),
      iter = 4,
      warmup = 2
    ),
    "`y` must be a Surv object, data frame, or matrix."
  )
})
