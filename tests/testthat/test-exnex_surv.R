testthat::local_edition(3)

expect_error_message <- function(expr, pattern) {
  err <- tryCatch(expr, error = function(e) e)
  testthat::expect_s3_class(err, "error")
  testthat::expect_match(conditionMessage(err), pattern, fixed = TRUE)
}

simulate_surv_data <- function(theta, sigma2, beta = NULL, n_per_group = 60, censor_min = 4, censor_max = 12) {
  groups <- rep(seq_along(theta), each = n_per_group)
  age_std <- rnorm(length(groups), mean = 0, sd = 1)
  mean_log_time <- rep(theta, each = n_per_group)

  if (!is.null(beta)) {
    mean_log_time <- mean_log_time + beta * age_std
  }

  log_time <- rnorm(length(groups), mean = mean_log_time, sd = sqrt(sigma2))
  true_time <- exp(log_time)
  censor_time <- runif(length(groups), min = censor_min, max = censor_max)
  time <- pmin(true_time, censor_time)
  event <- as.integer(true_time <= censor_time)

  out <- data.frame(
    time = time,
    event = event,
    group = factor(groups)
  )

  if (!is.null(beta)) {
    out$age_std <- age_std
  }

  out
}

trial_data <- data.frame(
  time = c(5, 8, 12, 9),
  event = c(1, 0, 1, 1),
  group = factor(c("A", "B", "A", "B")),
  age = c(60, 55, 62, 58)
)

test_that("exnex_surv works with formula and data.frame interfaces", {
  fit_formula <- exnex_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    priors = list(alpha = 1),
    iter = 6,
    warmup = 2,
    chains = 1,
    seed = 2719
  )

  fit_xy <- exnex_surv(
    x = trial_data[c("group", "age")],
    y = survival::Surv(trial_data$time, trial_data$event),
    priors = list(alpha = 1),
    iter = 6,
    warmup = 2,
    chains = 1,
    seed = 2719
  )

  expect_identical(
    class(fit_formula),
    c("exnex_surv", "hardhat_model", "hardhat_scalar")
  )
  expect_identical(
    class(fit_formula$blueprint),
    c("default_formula_blueprint", "formula_blueprint", "hardhat_blueprint")
  )
  expect_identical(
    class(fit_xy$blueprint),
    c("default_xy_blueprint", "xy_blueprint", "hardhat_blueprint")
  )
  expect_identical(dim(fit_formula$draws), c(4L, 4L))
  expect_identical(colnames(fit_formula$draws), c("theta_1", "theta_2", "beta_1", "sigma2"))
  expect_true(all(is.finite(as.matrix(fit_formula$draws))))
  expect_true(all(fit_formula$draws$sigma2 > 0))
  expect_identical(fit_formula$draws, fit_xy$draws)
  expect_identical(fit_formula$data$time, trial_data$time)
  expect_identical(fit_formula$data$event, trial_data$event)
  expect_identical(fit_formula$data$group, c(1, 2, 1, 2))
  expect_identical(fit_formula$data$X[, 1], trial_data$age)
  expect_identical(fit_formula$data$cov_names, "age")
  expect_equal(fit_formula$data$n, 4)
  expect_equal(fit_formula$data$n_groups, 2)
  expect_equal(fit_formula$data$n_covariates, 1)
})

test_that("exnex_surv handles models without covariates", {
  fit <- exnex_surv(
    survival::Surv(time, event) ~ group,
    data = trial_data,
    priors = list(),
    iter = 6,
    warmup = 2,
    chains = 1,
    seed = 5831
  )

  expect_identical(dim(fit$draws), c(4L, 3L))
  expect_identical(colnames(fit$draws), c("theta_1", "theta_2", "sigma2"))
  expect_true(all(is.finite(as.matrix(fit$draws))))
  expect_true(all(fit$draws$sigma2 > 0))
  expect_equal(fit$data$n_covariates, 0)
  expect_identical(fit$data$cov_names, character(0))
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
    iter = 6,
    warmup = 2,
    chains = 1,
    seed = 2719
  )

  expect_identical(colnames(fit$draws), c("theta_1", "theta_2", "beta_1", "beta_2", "sigma2"))
  expect_identical(colnames(fit$data$X), c("group_size", "age"))
  expect_identical(fit$data$cov_names, c("group_size", "age"))
  expect_identical(unname(fit$data$X[, "group_size"]), prefixed_data$group_size)
  expect_identical(unname(fit$data$X[, "age"]), prefixed_data$age)
  expect_equal(fit$data$n_covariates, 2)
})

test_that("exnex_surv roughly recovers simulated parameters with covariates", {
  set.seed(6841)
  sim_data <- simulate_surv_data(
    theta = c(1.2, 1.7, 2.1),
    sigma2 = 0.20,
    beta = -0.35,
    n_per_group = 60,
    censor_min = 4,
    censor_max = 12
  )

  fit <- exnex_surv(
    survival::Surv(time, event) ~ group + age_std,
    data = sim_data,
    iter = 500,
    warmup = 250,
    chains = 1,
    seed = 6841
  )

  draw_means <- colMeans(fit$draws)

  expect_true(mean(sim_data$event) > 0.5)
  expect_true(mean(sim_data$event) < 0.8)
  expect_equal(unname(draw_means[c("theta_1", "theta_2", "theta_3")]), c(1.2, 1.7, 2.1), tolerance = 0.35)
  expect_equal(unname(draw_means["beta_1"]), -0.35, tolerance = 0.20)
  expect_equal(unname(draw_means["sigma2"]), 0.20, tolerance = 0.20)
})

test_that("exnex_surv roughly recovers simulated parameters without covariates", {
  set.seed(5217)
  sim_data <- simulate_surv_data(
    theta = c(1.0, 1.8),
    sigma2 = 0.25,
    n_per_group = 70,
    censor_min = 3.5,
    censor_max = 10
  )

  fit <- exnex_surv(
    survival::Surv(time, event) ~ group,
    data = sim_data,
    iter = 500,
    warmup = 250,
    chains = 1,
    seed = 5217
  )

  draw_means <- colMeans(fit$draws)

  expect_true(mean(sim_data$event) > 0.6)
  expect_true(mean(sim_data$event) < 0.8)
  expect_equal(unname(draw_means[c("theta_1", "theta_2")]), c(1.0, 1.8), tolerance = 0.30)
  expect_equal(unname(draw_means["sigma2"]), 0.25, tolerance = 0.20)
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
