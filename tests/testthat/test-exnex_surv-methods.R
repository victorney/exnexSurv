testthat::local_edition(3)

test_that("summary.exnex_surv returns posterior table", {
  trial_data <- data.frame(
    time = c(5, 8, 12, 9),
    event = c(1, 0, 1, 1),
    group = factor(c("A", "B", "A", "B")),
    age = c(60, 55, 62, 58)
  )

  fit <- exnex_surv(
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
    c("parameter", "mean", "sd", "q05", "q50", "q95")
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

  fit <- exnex_surv(
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

  fit <- exnex_surv(
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

  fit <- exnex_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    iter = 8,
    warmup = 3,
    chains = 2,
    parallel_chains = 2,
    seed = 2719
  )

  png(filename = tempfile(fileext = ".png"), width = 600, height = 500)
  on.exit(dev.off(), add = TRUE)

  expect_invisible(plot(fit, parameters = "theta_1", ask = FALSE))
})
