#' Simulate basket-trial log-normal survival data
#'
#' Generates a synthetic dataset from the log-normal AFT model used by
#' `exnexSurv`:
#' \deqn{\log T_i = \theta_{g[i]} + X_i^\top\beta + \varepsilon,\qquad
#'   \varepsilon\sim\mathcal N(0,\sigma^2).}
#' Healthy baskets draw their location from a common baseline plus random noise;
#' "outlier" baskets (specified via `outlier_baskets`) have their location
#' shifted by `resist_delta`. The function is useful for simulation studies and
#' teaching examples.
#'
#' @param n Number of patients per basket. A scalar is replicated across all
#'   `K` baskets; a numeric vector of length `K` assigns a size to each basket
#'   individually.
#' @param K Number of baskets (default `9`).
#' @param beta Numeric vector of regression coefficients for the covariates.
#'   Length determines the number of covariates.
#' @param sigma Residual standard deviation (default `1.2`).
#' @param outlier_baskets Optional integer vector of basket indices (1 to K)
#'   whose true location is shifted away from the healthy population. These are
#'   the baskets the EXNEX model is designed to detect.
#' @param resist_delta Additive shift applied to the `theta` of outlier
#'   baskets (default `-0.8`).
#' @param censoring_rate Approximate proportion of censoring after
#'   `censor_upper`; if `NULL`, no censoring is applied.
#' @param censor_upper Upper bound of the censoring-time uniform distribution.
#' @param theta Baseline location for the healthy (non-outlier) baskets; a
#'   scalar used as the centre around which the healthy basket locations vary.
#' @param seed Optional seed for reproducibility.
#'
#' @return A `data.frame` with columns `time`, `event`, `group` (a factor), and
#'   one covariate column (`x1`, ...) per entry in `beta`. True parameter values
#'   are stored as attributes `true_theta`, `true_beta`, and `true_sigma`.
#' @export
simulate_data <- function(
  n = 30,
  K = 9,
  beta = c(0.5, -0.2),
  sigma = 1.2,
  outlier_baskets = NULL,
  resist_delta = -0.8,
  censoring_rate = NULL,
  censor_upper = NULL,
  theta = 0,
  seed = NULL
) {
  checkmate::assert_int(K, lower = 2)
  if (length(n) == 1L) {
    checkmate::assert_int(n, lower = 1)
    sizes <- rep(as.integer(n), K)
  } else if (length(n) == K) {
    checkmate::assert_integerish(n, lower = 1, any.missing = FALSE, len = K)
    sizes <- as.integer(n)
  } else {
    stop("`n` must be a scalar (recycled to all baskets) or a vector of length K = ",
         K, ".", call. = FALSE)
  }
  checkmate::assert_numeric(beta, any.missing = FALSE)
  checkmate::assert_number(sigma, lower = 0, finite = TRUE)
  checkmate::assert_numeric(outlier_baskets, lower = 1, upper = K,
                            any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_number(resist_delta, finite = TRUE)
  checkmate::assert_number(censoring_rate, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_number(censor_upper, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_number(theta, finite = TRUE)
  checkmate::assert_int(seed, lower = 1, upper = 2147483647, null.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  n_total <- sum(sizes)
  group <- rep(seq_len(K), times = sizes)

  # healthy baskets: baseline plus small random variation
  healthy_theta <- theta + rnorm(K, 0, 0.15)
  true_theta <- healthy_theta
  outlier_idx <- if (is.null(outlier_baskets)) integer(0) else outlier_baskets
  true_theta[outlier_idx] <- healthy_theta[outlier_idx] + resist_delta

  P <- length(beta)
  X <- matrix(rnorm(n_total * P, mean = 0, sd = 1), nrow = n_total, ncol = P)

  eta <- true_theta[group]
  if (P > 0) eta <- eta + X %*% beta

  log_t <- rnorm(n_total, mean = eta, sd = sigma)
  true_time <- exp(log_t)

  event <- rep(1L, n_total)
  time <- true_time
  if (!is.null(censoring_rate)) {
    if (is.null(censor_upper)) censor_upper <- as.numeric(stats::quantile(true_time, 0.8))
    # adapt upper bound until approximately the requested censoring rate
    censor_time <- runif(n_total, 0, censor_upper)
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
