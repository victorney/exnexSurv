#' Simulate basket-trial log-normal survival data
#'
#' Generates a synthetic dataset from the log-normal AFT model used by
#' `exnexSurv`:
#' \deqn{\log T_i = \theta_{g[i]} + X_i^\top\beta + \varepsilon,\qquad
#'   \varepsilon\sim\mathcal N(0,\sigma^2).}
#' Healthy (non-resistant) baskets draw their location from a common baseline
#' plus random noise; "resistant" baskets (the ones the EXNEX model is designed
#' to detect) have their location shifted by `resist_delta`. The function is
#' useful for simulation studies and teaching examples.
#'
#' @param n_each Number of patients in each basket.
#' @param beta Numeric vector of regression coefficients for the covariates.
#'   Length determines the number of covariates.
#' @param sigma Residual standard deviation (default `1.2`).
#' @param resistant Optional integer vector of basket indices to make
#'   "resistant" (their true `theta` is shifted away from the population).
#' @param resist_delta Additive shift applied to the `theta` of resistant
#'   baskets (default `-0.8`).
#' @param censoring_rate Approximate proportion of censoring after
#'   `censor_upper`; if `NULL`, no censoring is applied.
#' @param censor_upper Upper bound of the censoring-time uniform distribution.
#' @param theta Null-posterior location for the non-resistant baskets; a scalar
#'   used as a base around which healthy baskets vary.
#' @param seed Optional seed for reproducibility.
#'
#' @return A `data.frame` with columns `time`, `event`, `group` (a factor), and
#'   one covariate column (`x1`, ...) per entry in `beta`. True parameter values
#'   are stored as attributes `true_theta`, `true_beta`, and `true_sigma`.
#' @export
simulate_data <- function(
  n_each = 30,
  beta = c(0.5, -0.2),
  sigma = 1.2,
  resistant = NULL,
  resist_delta = -0.8,
  censoring_rate = NULL,
  censor_upper = NULL,
  theta = 0,
  seed = NULL
) {
  checkmate::assert_int(n_each, lower = 1)
  checkmate::assert_numeric(beta, any.missing = FALSE)
  checkmate::assert_number(sigma, lower = 0, finite = TRUE)
  checkmate::assert_numeric(resistant, lower = 1, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_number(resist_delta, finite = TRUE)
  checkmate::assert_number(censoring_rate, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_number(censor_upper, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_number(theta, finite = TRUE)
  checkmate::assert_int(seed, lower = 1, upper = 2147483647, null.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  K <- 9L  # fixed number of baskets for this generator
  if (!is.null(resistant) && max(resistant) > K) {
    stop("`resistant` indices must be <= ", K, call. = FALSE)
  }

  # allocate per-basket sizes from 30 down to 5 (like the case study)
  sizes <- round(seq(30, 5, length.out = K))
  n <- sum(sizes)
  group <- rep(seq_len(K), times = sizes)

  # healthy baskets: baseline plus small random variation
  healthy_theta <- theta + rnorm(K, 0, 0.15)
  true_theta <- healthy_theta
  resist <- if (is.null(resistant)) integer(0) else resistant
  true_theta[resist] <- healthy_theta[resist] + resist_delta

  P <- length(beta)
  X <- matrix(rnorm(n * P, mean = 0, sd = 1), nrow = n, ncol = P)

  eta <- true_theta[group]
  if (P > 0) eta <- eta + X %*% beta

  log_t <- rnorm(n, mean = eta, sd = sigma)
  true_time <- exp(log_t)

  event <- rep(1L, n)
  time <- true_time
  if (!is.null(censoring_rate)) {
    if (is.null(censor_upper)) censor_upper <- as.numeric(stats::quantile(true_time, 0.8))
    # adapt upper bound until approximately the requested censoring rate
    censor_time <- runif(n, 0, censor_upper)
    event <- as.integer(true_time <= censor_time)
    time <- pmin(true_time, censor_time)
  }

  cov_df <- if (P > 0) {
    stats::setNames(as.data.frame(X), paste0("x", seq_len(P)))
  } else {
    NULL
  }

  out <- data.frame(time = time, event = event, group = factor(group))
  if (P > 0) out <- cbind(out, cov_df)

  attr(out, "true_theta") <- true_theta
  attr(out, "true_beta") <- beta
  attr(out, "true_sigma") <- sigma

  out
}
