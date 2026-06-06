testthat::local_edition(3)

trial_data_local <- data.frame(
  time = c(5, 8, 12, 9),
  event = c(1, 0, 1, 1),
  group = factor(c("A", "B", "A", "B")),
  age = c(60, 55, 62, 58)
)

expect_error_message <- function(expr, pattern) {
  err <- tryCatch(expr, error = function(e) e)
  testthat::expect_s3_class(err, "error")
  testthat::expect_match(conditionMessage(err), pattern, fixed = TRUE)
}

make_surv <- function() {
  survival::Surv(c(5, 8, 12), c(1, 0, 1))
}

test_that("bridge helpers extract survival vectors", {
  surv_obj <- make_surv()
  outcomes_df <- data.frame(time = c(5, 8, 12), event = c(1, 0, 1))

  expect_identical(
    exnexSurv:::.extract_time_vector(list(surv_obj)),
    c(5, 8, 12)
  )
  expect_identical(
    exnexSurv:::.extract_event_vector(list(surv_obj)),
    c(1, 0, 1)
  )
  expect_identical(exnexSurv:::.extract_time_vector(outcomes_df), c(5, 8, 12))
  expect_identical(exnexSurv:::.extract_event_vector(outcomes_df), c(1, 0, 1))
})

test_that("bridge helpers validate malformed outcomes", {
  expect_error_message(
    exnexSurv:::.extract_time_vector(list(data.frame(time = 1))),
    "Outcome must contain time and event status."
  )

  expect_error_message(
    exnexSurv:::.extract_event_vector(list(data.frame(time = 1))),
    "Outcome must contain time and event status."
  )
})

test_that("new_exnex_surv enforces structural invariants", {
  data <- list(
    time = c(5, 8, 12),
    event = c(1, 0, 1),
    group = c(1, 2, 1),
    X = matrix(c(60, 55, 62), ncol = 1),
    n = 3,
    n_groups = 2,
    n_covariates = 1,
    cov_names = "age"
  )
  draws <- data.frame(
    theta_1 = c(0, 0),
    theta_2 = c(0, 0),
    beta_1 = c(0, 0),
    sigma = c(1, 1)
  )

  fit <- exnex_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data_local,
    priors = list(alpha = 1),
    iter = 4,
    warmup = 2,
    chains = 1
  )

  rebuilt <- exnexSurv:::new_exnex_surv(
    draws = draws,
    data = data,
    priors = list(alpha = 1),
    iter = 4,
    warmup = 2,
    chains = 1,
    blueprint = fit$blueprint
  )

  expect_s3_class(rebuilt, "exnex_surv")
  expect_identical(rebuilt$draws, draws)
  expect_identical(rebuilt$data$cov_names, "age")

  expect_error_message(
    exnexSurv:::new_exnex_surv(
      draws = draws[1, , drop = FALSE],
      data = data,
      priors = list(),
      iter = 4,
      warmup = 2,
      chains = 1,
      blueprint = fit$blueprint
    ),
    "Number of draw rows"
  )

  expect_error_message(
    exnexSurv:::new_exnex_surv(
      draws = draws[, 1:3],
      data = data,
      priors = list(),
      iter = 4,
      warmup = 2,
      chains = 1,
      blueprint = fit$blueprint
    ),
    "Number of columns in draws"
  )
})

test_that("public wrappers keep rejecting invalid structural inputs", {
  expect_error_message(
    exnex_surv(
      survival::Surv(time, event) ~ group,
      data = transform(trial_data_local, time = c(5, -8, 12, 9)),
      priors = list(),
      iter = 4,
      warmup = 2
    ),
    "All survival times must be positive."
  )

  expect_error(
    suppressWarnings(
      exnex_surv(
        survival::Surv(time, event) ~ group,
        data = transform(trial_data_local, event = c(1, 0, 2, 1)),
        priors = list(),
        iter = 4,
        warmup = 2
      )
    )
  )

  expect_error_message(
    exnex_surv(
      x = trial_data_local[c("group", "age")],
      y = survival::Surv(trial_data_local$time, trial_data_local$event),
      priors = list(),
      iter = 4,
      warmup = 2,
      group_col = "missing"
    ),
    "Column 'missing' not found in predictors."
  )
})
