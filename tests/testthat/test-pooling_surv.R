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
  fit_formula <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    priors = list(alpha = 1),
    iter = 6,
    warmup = 2,
    chains = 1,
    seed = 2719
  )

  fit_xy <- pooling_surv(
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
    c("pooling_surv", "hardhat_model", "hardhat_scalar")
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
  fit <- pooling_surv(
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

  fit <- pooling_surv(
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

  fit <- pooling_surv(
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

  fit <- pooling_surv(
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
    pooling_surv(
      survival::Surv(time, event) ~ group,
      data = trial_data,
      priors = list(),
      iter = 4,
      warmup = 4
    ),
    "`warmup` (4) must be less than `iter` (4)."
  )

  fit_multichain <- pooling_surv(
    survival::Surv(time, event) ~ group,
    data = trial_data,
    priors = list(),
    iter = 4,
    warmup = 2,
    chains = 2,
    seed = 2719
  )

  expect_identical(dim(fit_multichain$draws), c(4L, 3L)) 

  expect_error_message(
    pooling_surv(
      survival::Surv(time, event) ~ group,
      data = trial_data,
      priors = list(),
      iter = 4,
      warmup = 2,
      chains = 2,
      parallel_chains = "yes"
    ),
    "Must be of type 'logical flag'"
  )

  expect_error_message(
    pooling_surv(
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
    pooling_surv(
      x = trial_data[c("group", "age")],
      y = trial_data$time,
      priors = list(),
      iter = 4,
      warmup = 2
    ),
    "`y` must be a Surv object, data frame, or matrix."
  )
})

test_that("exnex_surv supports parallel chain execution", {
  fit_seq <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    priors = list(alpha = 1),
    iter = 8,
    warmup = 3,
    chains = 2,
    parallel_chains = FALSE,
    seed = 2719,
    verbose = FALSE
  )

  fit_parallel <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    priors = list(alpha = 1),
    iter = 8,
    warmup = 3,
    chains = 2,
    parallel_chains = TRUE,
    seed = 2719,
    verbose = FALSE
  )

  expect_identical(dim(fit_seq$draws), c(10L, 4L))
  expect_identical(fit_seq$draws, fit_parallel$draws)
  expect_identical(fit_seq$data$chain_seeds, fit_parallel$data$chain_seeds)
})

test_that("pooling_surv validates the pooling argument", {
  expect_error_message(
    pooling_surv(
      survival::Surv(time, event) ~ group,
      data = trial_data,
      pooling = "banana",
      iter = 4,
      warmup = 2
    ),
    "Must be element of set"
  )
})

test_that("pooling = 'complete' shares one basket effect across groups", {
  set.seed(2719)
  d <- simulate_surv_data(theta = c(1.5, 1.7), sigma2 = 0.5, n_per_group = 40)

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group,
    data = d,
    pooling = "complete",
    iter = 200,
    warmup = 100,
    chains = 1,
    seed = 2719
  )

  expect_s3_class(fit, "pooling_surv")
  expect_identical(fit$pooling, "complete")
  expect_identical(dim(fit$draws), c(100L, 3L))
  expect_identical(colnames(fit$draws), c("theta_1", "theta_2", "sigma2"))

  # All theta columns hold identical draws from the shared effect
  expect_identical(fit$draws$theta_1, fit$draws$theta_2)
})

test_that("complete pooling matches the analytic conjugate posterior (no censoring)", {
  set.seed(42)
  theta <- c(1.2, 1.6)
  sigma2 <- 0.4
  n_per_group <- 50
  log_time <- c(
    rnorm(n_per_group, theta[1], sqrt(sigma2)),
    rnorm(n_per_group, theta[2], sqrt(sigma2))
  )
  d <- data.frame(
    time = exp(log_time),
    event = rep(1, 2 * n_per_group),
    group = factor(rep(1:2, each = n_per_group))
  )

  # Prior N(m_mu = 0, v_mu = 1e4); analytic marginal posterior of theta
  y_bar <- mean(log_time)
  prior_var <- 1e4
  post_var <- 1 / (1 / prior_var + (2 * n_per_group) / sigma2)
  post_mean <- post_var * (y_bar * (2 * n_per_group) / sigma2)

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group,
    data = d,
    pooling = "complete",
    iter = 2000,
    warmup = 400,
    chains = 1,
    seed = 11
  )

  expect_equal(mean(fit$draws$theta_1), post_mean, tolerance = 0.02)
})

test_that("no pooling matches analytic posteriors for each basket", {
  set.seed(42)
  theta <- c(1.0, 2.0)
  sigma2 <- 0.5
  n_per_group <- 60
  log_time <- c(
    rnorm(n_per_group, theta[1], sqrt(sigma2)),
    rnorm(n_per_group, theta[2], sqrt(sigma2))
  )
  d <- data.frame(
    time = exp(log_time),
    event = rep(1, 2 * n_per_group),
    group = factor(rep(1:2, each = n_per_group))
  )

  # No pooling: theta_j ~ N(m_nex_j = 0, v_nex_j = 1e4), independent
  post_var <- 1 / (1 / 1e4 + n_per_group / sigma2)
  post_means <- c(
    mean(log_time[1:n_per_group]) * (n_per_group / sigma2),
    mean(log_time[-(1:n_per_group)]) * (n_per_group / sigma2)
  ) * post_var

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group,
    data = d,
    pooling = "none",
    iter = 2000,
    warmup = 400,
    chains = 1,
    seed = 11
  )

  expect_equal(mean(fit$draws$theta_1), post_means[1], tolerance = 0.02)
  expect_equal(mean(fit$draws$theta_2), post_means[2], tolerance = 0.02)
})

test_that("exnex mode is unchanged by default pooling", {
  set.seed(2719)
  d <- simulate_surv_data(theta = c(1.1, 1.9), sigma2 = 0.4, n_per_group = 30)

  fit_default <- pooling_surv(
    survival::Surv(time, event) ~ group,
    data = d,
    iter = 60,
    warmup = 20,
    chains = 1,
    seed = 2719
  )
  fit_exnex <- pooling_surv(
    survival::Surv(time, event) ~ group,
    data = d,
    pooling = "exnex",
    iter = 60,
    warmup = 20,
    chains = 1,
    seed = 2719
  )

  expect_identical(fit_default$draws, fit_exnex$draws)
})

test_that("priors unused in the selected pooling mode trigger one warning", {
  set.seed(2719)
  d <- simulate_surv_data(theta = c(1.5, 1.7), sigma2 = 0.5, n_per_group = 20)

  expect_warning(
    fit_complete <- pooling_surv(
      survival::Surv(time, event) ~ group,
      data = d,
      pooling = "complete",
      priors = list(p_mix = 0.7, a_tau = 3, b_tau = 4, v_beta = 1e4),
      iter = 20,
      warmup = 5,
      chains = 1,
      seed = 2719
    ),
    "are ignored with pooling = \"complete\"",
    fixed = TRUE
  )

  expect_warning(
    pooling_surv(
      survival::Surv(time, event) ~ group,
      data = d,
      pooling = "none",
      priors = list(m_mu = 1, v_mu = 4, p_mix = 0.5),
      iter = 20,
      warmup = 5,
      chains = 1,
      seed = 2719
    ),
    "are ignored with pooling = \"none\"",
    fixed = TRUE
  )

  # Relevant priors alone raise no warning
  expect_no_warning(
    pooling_surv(
      survival::Surv(time, event) ~ group,
      data = d,
      pooling = "none",
      priors = list(m_nex = 1, v_nex = 4),
      iter = 20,
      warmup = 5,
      chains = 1,
      seed = 2719
    )
  )

  expect_s3_class(fit_complete, "pooling_surv")
})

test_that("exnex_surv() alias forces pooling = exnex with a deprecation message", {
  set.seed(2719)

  expect_message(
    fit_alias <- exnex_surv(
      survival::Surv(time, event) ~ group + age,
      data = trial_data,
      priors = list(alpha = 1),
      iter = 6,
      warmup = 2,
      chains = 1,
      seed = 2719
    ),
    "exnex_surv() is deprecated; use pooling_surv() instead",
    fixed = TRUE
  )

  fit_new <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    priors = list(alpha = 1),
    pooling = "exnex",
    iter = 6,
    warmup = 2,
    chains = 1,
    seed = 2719
  )

  expect_identical(fit_alias$draws, fit_new$draws)
  expect_identical(fit_alias$pooling, "exnex")
})

test_that("verbose shows progress and can be silenced", {
  # verbose = TRUE draws a progress bar (progressr) and must not error
  fit_v <- suppressMessages(pooling_surv(
    survival::Surv(time, event) ~ group,
    data = trial_data,
    iter = 30,
    warmup = 10,
    chains = 1,
    seed = 2719
  ))
  expect_s3_class(fit_v, "pooling_surv")

  expect_silent(
    fit_nv <- pooling_surv(
      survival::Surv(time, event) ~ group,
      data = trial_data,
      iter = 30,
      warmup = 10,
      chains = 1,
      seed = 2719,
      verbose = FALSE
    )
  )
  expect_s3_class(fit_nv, "pooling_surv")
})

test_that("pooling_surv thins draws and stores the interval", {
  fit_full <- pooling_surv(
    survival::Surv(time, event) ~ group,
    data = trial_data,
    iter = 30,
    warmup = 10,
    thin = 1,
    chains = 2,
    seed = 2719,
    verbose = FALSE
  )
  fit_thin <- pooling_surv(
    survival::Surv(time, event) ~ group,
    data = trial_data,
    iter = 30,
    warmup = 10,
    thin = 2,
    chains = 2,
    seed = 2719,
    verbose = FALSE
  )

  expect_identical(fit_thin$thin, 2L)
  expect_identical(nrow(fit_thin$draws), 20L)   # ceiling(20 / 2) * 2 chains
  expect_identical(nrow(fit_full$draws), 40L)   # 20 * 2 chains

  # thinning keeps every other post-warmup draw of the same chain
  chain1_full <- fit_full$draws[1:20, , drop = FALSE]
  chain1_thin <- fit_thin$draws[1:10, , drop = FALSE]
  expect_equal(
    unname(as.matrix(chain1_thin)),
    unname(as.matrix(chain1_full[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19), , drop = FALSE]))
  )

  fit_thin_par <- pooling_surv(
    survival::Surv(time, event) ~ group,
    data = trial_data,
    iter = 30,
    warmup = 10,
    thin = 2,
    chains = 2,
    parallel_chains = TRUE,
    seed = 2719,
    verbose = FALSE
  )
  expect_equal(fit_thin_par$draws, fit_thin$draws)
  expect_identical(fit_thin_par$thin, 2L)

  expect_error_message(
    pooling_surv(
      survival::Surv(time, event) ~ group,
      data = trial_data,
      iter = 30,
      warmup = 10,
      thin = 0,
      seed = 2719,
      verbose = FALSE
    ),
    "Element 1 is not >= 1"
  )
})
