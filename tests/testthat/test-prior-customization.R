expect_error_message <- function(expr, pattern) {
  err <- tryCatch(expr, error = function(e) e)
  testthat::expect_s3_class(err, "error")
  testthat::expect_match(conditionMessage(err), pattern, fixed = TRUE)
}

base_args <- function(priors = list(), verbose = FALSE,
    iter = 20, warmup = 10) {
  list(
    time = c(5, 8, 12, 9, 7, 11),
    event = c(1, 0, 1, 1, 0, 1),
    group = c(1, 2, 1, 2, 1, 2),
    X = matrix(c(60, 55, 62, 58, 61, 59), ncol = 1),
    priors = priors,
    pooling = "exnex",
    verbose = FALSE,
    iter = iter,
    warmup = warmup,
    chains = 1,
    chain_label = ""
  )
}

test_that("absent and unknown priors keep defaults and reproducibility", {
  set.seed(101)
  r_default <- do.call(exnexSurv:::cpp_exnex_gibbs, base_args(list()))
  set.seed(101)
  r_unknown <- do.call(
    exnexSurv:::cpp_exnex_gibbs,
    base_args(list(alpha = 1))
  )

  expect_identical(r_default$draws, r_unknown$draws)
  expect_identical(r_unknown$priors, list(alpha = 1))

  resolved <- r_default$diagnostics$resolved_priors
  expect_equal(resolved$a_sigma, 2.0)
  expect_equal(resolved$b_sigma, 2.0)
  expect_equal(resolved$a_tau, 2.0)
  expect_equal(resolved$b_tau, 2.0)
  expect_equal(resolved$p_mix, 0.5)
  expect_equal(resolved$m_mu, 0.0)
  expect_equal(resolved$v_mu, 1e4)
  expect_equal(resolved$m_nex, 0.0)
  expect_equal(resolved$v_nex, 1e4)
  expect_equal(resolved$v_beta, 1e4)
})

test_that("custom priors are read and reported", {
  set.seed(202)
  res <- do.call(
    exnexSurv:::cpp_exnex_gibbs,
    base_args(list(
      a_sigma = 3, b_sigma = 3, a_tau = 3, b_tau = 3,
      p_mix = 0.7, m_mu = 0.5, v_mu = 10, m_nex = 1, v_nex = 5
    ))
  )

  resolved <- res$diagnostics$resolved_priors
  expect_equal(resolved$a_sigma, 3)
  expect_equal(resolved$b_sigma, 3)
  expect_equal(resolved$a_tau, 3)
  expect_equal(resolved$b_tau, 3)
  expect_equal(resolved$p_mix, 0.7)
  expect_equal(resolved$m_mu, 0.5)
  expect_equal(resolved$v_mu, 10)
  expect_equal(resolved$m_nex, 1)
  expect_equal(resolved$v_nex, 5)
  expect_equal(resolved$v_beta, 1e4)
  expect_true(all(is.finite(res$draws)))
})

test_that("custom priors are reproducible", {
  args <- base_args(list(a_tau = 3, p_mix = 0.7, v_nex = 10))
  set.seed(303)
  r1 <- do.call(exnexSurv:::cpp_exnex_gibbs, args)
  set.seed(303)
  r2 <- do.call(exnexSurv:::cpp_exnex_gibbs, args)
  expect_identical(r1$draws, r2$draws)
})

test_that("invalid priors raise clear errors", {
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(a_sigma = 0))
    ),
    "Inverse-gamma prior parameters must be positive."
  )
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(b_tau = -1))
    ),
    "Inverse-gamma prior parameters must be positive."
  )
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(p_mix = 1.5))
    ),
    "must lie strictly between 0 and 1"
  )
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(p_mix = 0))
    ),
    "must lie strictly between 0 and 1"
  )
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(v_mu = NaN))
    ),
    "must be finite"
  )
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(m_mu = "a"))
    ),
    "must be a single numeric value"
  )
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(a_tau = c(1, 2)))
    ),
    "must be a single numeric value"
  )
})

test_that("mixture weight near one keeps the sampler stable", {
  set.seed(404)
  res <- do.call(
    exnexSurv:::cpp_exnex_gibbs,
    base_args(list(p_mix = 0.9999999))
  )
  expect_equal(res$diagnostics$resolved_priors$p_mix, 0.9999999)
  expect_true(all(is.finite(res$draws)))
})

test_that("basket-specific vector priors are accepted and reported", {
  set.seed(606)
  res <- do.call(
    exnexSurv:::cpp_exnex_gibbs,
    base_args(list(
      p_mix = c(0.6, 0.4),
      m_nex = c(-1, 2),
      v_nex = c(1, 100)
    ))
  )

  resolved <- res$diagnostics$resolved_priors
  expect_equal(resolved$p_mix, c(0.6, 0.4))
  expect_equal(resolved$m_nex, c(-1, 2))
  expect_equal(resolved$v_nex, c(1, 100))
  expect_true(all(is.finite(res$draws)))
})

test_that("vector priors pull basket effects toward their own NEX component", {
  set.seed(707)
  # Basket 1 and 3 are exchangeable around log-time 5; basket 2 has data
  # around 0 but is forced into its own tight NEX prior centered at 0.
  time <- exp(c(rep(5, 20), rep(0, 20), rep(5, 20)) + rnorm(60, 0, 0.2))
  event <- rep(1, 60)
  group <- rep(1:3, each = 20)
  res <- exnexSurv:::cpp_exnex_gibbs(
    time = time,
    event = event,
    group = group,
    X = matrix(nrow = 60, ncol = 0),
    priors = list(
      p_mix = c(0.9999999, 1e-6, 0.9999999),
      m_nex = c(0, 0, 0),
      v_nex = c(1e4, 1e-4, 1e4)
    ),
    pooling = "exnex",
    verbose = FALSE,
    iter = 1000,
    warmup = 500,
    chains = 1,
    chain_label = ""
  )
  posterior_means <- colMeans(res$draws[, 1:3])
  expect_true(abs(posterior_means[2]) < 0.3)      # pulled to 0 by NEX prior
  expect_true(abs(posterior_means[1] - 5) < 0.5)  # exchangeable basket 1
  expect_true(abs(posterior_means[3] - 5) < 0.5)  # exchangeable basket 3
})

test_that("vector priors of the wrong length raise clear errors", {
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(p_mix = c(0.5, 0.5, 0.5)))
    ),
    "must be a single numeric value or a numeric vector of length K"
  )
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(v_nex = c(1, 2, 3)))
    ),
    "must be a single numeric value or a numeric vector of length K"
  )
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(p_mix = c(0.5, 1.5)))
    ),
    "must lie strictly between 0 and 1 for every basket"
  )
  expect_error_message(
    do.call(
      exnexSurv:::cpp_exnex_gibbs,
      base_args(list(v_nex = c(1, -1)))
    ),
    "must be positive for every basket"
  )
})

test_that("v_beta is configurable and reported", {
  set.seed(808)
  res <- do.call(
    exnexSurv:::cpp_exnex_gibbs,
    base_args(list(v_beta = 0.5))
  )
  expect_equal(res$diagnostics$resolved_priors$v_beta, 0.5)
  expect_true(all(is.finite(res$draws)))

  expect_error_message(
    do.call(exnexSurv:::cpp_exnex_gibbs, base_args(list(v_beta = 0))),
    "Prior variances must be positive."
  )
})

test_that("informative nonexchangeable prior pulls basket effects toward its mean", {
  set.seed(505)
  time <- exp(5 + rnorm(60, 0, 0.2))
  event <- rep(1, 60)
  group <- rep(1:3, each = 20)
  res <- exnexSurv:::cpp_exnex_gibbs(
    time = time,
    event = event,
    group = group,
    X = matrix(nrow = 60, ncol = 0),
    priors = list(p_mix = 1e-6, m_nex = 5, v_nex = 1e-4),
    pooling = "exnex",
    verbose = FALSE,
    iter = 1000,
    warmup = 500,
    chains = 1,
    chain_label = ""
  )
  posterior_means <- colMeans(res$draws[, 1:3])
  expect_true(all(abs(posterior_means - 5) < 0.3))
})

test_that("pooling_surv() exposes resolved priors", {
  trial_data <- data.frame(
    time = c(5, 8, 12, 9, 7, 11),
    event = c(1, 0, 1, 1, 0, 1),
    group = c(1, 2, 1, 2, 1, 2),
    age = c(60, 55, 62, 58, 61, 59)
  )
  fit <- pooling_surv(
    survival::Surv(time, event) ~ group + age,
    data = trial_data,
    priors = list(p_mix = 0.7, v_nex = 10, a_tau = 3, b_tau = 3),
    verbose = FALSE,
    iter = 20,
    warmup = 10,
    chains = 1,
    chain_label = ""
  )

  resolved <- fit$resolved_priors
  expect_equal(resolved$p_mix, 0.7)
  expect_equal(resolved$v_nex, 10)
  expect_equal(resolved$a_tau, 3)
  expect_equal(resolved$b_tau, 3)
  expect_equal(resolved$a_sigma, 2.0)
  expect_equal(resolved$m_nex, 0.0)
  expect_equal(resolved$v_beta, 1e4)
})
