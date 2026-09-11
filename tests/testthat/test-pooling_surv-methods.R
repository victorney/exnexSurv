testthat::local_edition(3)

test_that("summary.exnex_surv returns posterior table", {
  trial_data <- data.frame(
    time = c(5, 8, 12, 9),
    event = c(1, 0, 1, 1),
    group = factor(c("A", "B", "A", "B")),
    age = c(60, 55, 62, 58)
  )

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    iter = 8,
    warmup = 3,
    chains = 1,
    seed = 2719
  )

  out <- summary(fit)

  expect_s3_class(out, "data.frame")
  expect_identical(
    colnames(out),
    c("parameter", "mean", "sd", "q05", "q50", "q95", "rhat", "ess_bulk", "ess_tail")
  )
  expect_identical(out$parameter, colnames(fit$draws))
})

test_that("print.exnex_surv works without plotting", {
  trial_data <- data.frame(
    time = c(5, 8, 12, 9),
    event = c(1, 0, 1, 1),
    group = factor(c("A", "B", "A", "B")),
    age = c(60, 55, 62, 58)
  )

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    iter = 8,
    warmup = 3,
    chains = 1,
    seed = 3119
  )

  expect_invisible(print(fit, show_trace = FALSE))
})

test_that("plot.exnex_surv creates per-parameter traceplot", {
  testthat::skip_if_not_installed("bayesplot")
  testthat::skip_if_not_installed("viridisLite")

  trial_data <- data.frame(
    time = c(5, 8, 12, 9),
    event = c(1, 0, 1, 1),
    group = factor(c("A", "B", "A", "B")),
    age = c(60, 55, 62, 58)
  )

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    iter = 10,
    warmup = 4,
    chains = 1,
    seed = 1447
  )

  png(filename = tempfile(fileext = ".png"), width = 600, height = 500)
  on.exit(dev.off(), add = TRUE)

  expect_invisible(plot(fit, parameters = "theta_1", ask = FALSE))
})

test_that("plot.exnex_surv overlays multiple chains per parameter", {
  testthat::skip_if_not_installed("bayesplot")
  testthat::skip_if_not_installed("viridisLite")

  trial_data <- data.frame(
    time = c(5, 8, 12, 9),
    event = c(1, 0, 1, 1),
    group = factor(c("A", "B", "A", "B")),
    age = c(60, 55, 62, 58)
  )

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    iter = 8,
    warmup = 3,
    chains = 2,
    parallel_chains = TRUE,
    seed = 2719
  )

  png(filename = tempfile(fileext = ".png"), width = 600, height = 500)
  on.exit(dev.off(), add = TRUE)

  expect_invisible(plot(fit, parameters = "theta_1", ask = FALSE))
})

test_that("convergence diagnostics summary include rhat and ess", {
  trial_data <- data.frame(
    time = c(5, 8, 12, 9),
    event = c(1, 0, 1, 1),
    group = factor(c("A", "B", "A", "B")),
    age = c(60, 55, 62, 58)
  )

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    iter = 8,
    warmup = 3,
    chains = 2,
    parallel_chains = TRUE,
    seed = 2719,
    verbose = FALSE
  )

  cd <- convergence_diagnostics(fit)
  expect_s3_class(cd, "data.frame")
  expect_identical(
    colnames(cd),
    c("parameter", "rhat", "ess_bulk", "ess_tail")
  )
  expect_identical(cd$parameter, colnames(fit$draws))
  expect_true(all(is.finite(cd$rhat)))

  summ <- summary(fit)
  expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% colnames(summ)))
})

test_that("rank_probabilities sums to one and ranks are dense", {
  trial_data <- data.frame(
    time = c(5, 8, 12, 9),
    event = c(1, 0, 1, 1),
    group = factor(c("A", "B", "A", "B"))
  )

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group,
    data = trial_data,
    iter = 8,
    warmup = 3,
    chains = 1,
    seed = 2719,
    verbose = FALSE
  )

  rp <- rank_probabilities(fit)
  expect_s3_class(rp, "data.frame")
  expect_identical(colnames(rp), c("group", "rank", "prob"))
  expect_equal(sum(rp$prob[rp$rank == 1]), 1)
  expect_setequal(unique(rp$rank), seq_len(2))
  expect_true(all(rp$prob >= 0 & rp$prob <= 1))
})

test_that("posterior_predictive_check and shrinkage_plot run silently", {
  trial_data <- data.frame(
    time = c(5, 8, 12, 9),
    event = c(1, 0, 1, 1),
    group = factor(c("A", "B", "A", "B")),
    age = c(60, 55, 62, 58)
  )

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    iter = 8,
    warmup = 3,
    chains = 1,
    seed = 2719,
    verbose = FALSE
  )

  suppressWarnings({pp <- posterior_predictive_check(fit, plot = FALSE, n_replicates = 5)})
  expect_equal(
    names(pp),
    c("posterior", "observed")
  )
  expect_true(all(c("group", "time", "surv") %in% colnames(pp$observed)))
  expect_true(all(c("group", "time", "mean") %in% colnames(pp$posterior)))

  sp <- suppressWarnings(shrinkage_plot(fit))
  expect_true(all(
    c("group", "posterior_mean", "posterior_lwr", "posterior_upr", "unpooled") %in% colnames(sp)
  ))
})

test_that("convergence diagnostics match manual rank-normalized split formulas", {
  skip_if_not_installed("posterior")

  trial_data <- data.frame(
    time = c(5, 8, 12, 9),
    event = c(1, 0, 1, 1),
    group = factor(c("A", "B", "A", "B")),
    age = c(60, 55, 62, 58)
  )

  fit <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    iter = 30,
    warmup = 10,
    chains = 2,
    seed = 2719
  )

  cd <- convergence_diagnostics(fit)
  arr <- exnexSurv:::.draws_array_pooling(fit)
  da <- posterior::as_draws_array(arr)

  expected_rhat <- c()
  for (v in dimnames(arr)[[3]]) {
    sub <- posterior::subset_draws(da, variable = v)
    expected_rhat <- c(expected_rhat, as.numeric(posterior::rhat(sub)))
  }
  expect_identical(cd$rhat, expected_rhat)
  expect_true(all(abs(cd$ess_bulk) > 0))
})
