#' Pointwise log-likelihood and WAIC from an exnex_surv fit
#'
#' Computes the Watanabe-Akaike Information Criterion (WAIC) for a fitted
#' `exnex_surv` model. WAIC is a fully Bayesian information criterion that uses
#' the posterior draws of the log-likelihood and adds a penalty for effective
#' number of parameters \eqn{p_{waic}}.
#'
#' For each observation \eqn{i}, with posterior draws indexed by
#' \eqn{s=1,\ldots,S}, we compute the pointwise log-likelihood
#' \eqn{\log p(y_i \mid \theta^{(s)})}, where the log-normal AFT model gives
#' for an observed event (\eqn{\delta_i=1})
#' \deqn{\log p(y_i \mid \theta^{(s)}) =
#'   \log\phi\!\left(\frac{\log t_i-\eta_i^{(s)}}{\sigma^{(s)}}\right)
#'   - \log(\sigma^{(s)} t_i),}
#' and for a censored observation (\eqn{\delta_i=0})
#' \deqn{\log p(y_i \mid \theta^{(s)}) =
#'   \log\!\left[1-\Phi\!\left(\frac{\log t_i-\eta_i^{(s)}}{\sigma^{(s)}}\right)\right],}
#' where \eqn{\eta_i^{(s)}} is the linear predictor built from the draws.
#'
#' WAIC is then
#' \deqn{lpd = \sum_i \log\!\left(\frac{1}{S}\sum_s \exp\log p(y_i\mid\theta^{(s)})\right),}
#' \deqn{p_{waic} = \sum_i \mathrm{var}_s\big(\log p(y_i\mid\theta^{(s)})\big),}
#' \deqn{\mathrm{WAIC} = -2(lpd - p_{waic}).}
#' Rows with extreme leverage (`p_waic` above a large threshold) may flag issues;
#' the function returns a `pointwise` matrix so users can investigate.
#'
#' @param fit A fitted `exnex_surv` object.
#' @param ... Unused.
#'
#' @return A named list with elements `waic`, `se_waic`, `lpd`, `p_waic`,
#'   `elpd_waic`, and `pointwise` (an `n` by 3 matrix with columns `lpd`,
#'   `p_waic`, `waic`).
#' @export
compute_waic <- function(fit, ...) {
  checkmate::assert_class(fit, "exnex_surv")

  draws <- fit$draws
  data <- fit$data
  K <- data$n_groups
  P <- data$n_covariates

  theta_matrix <- as.matrix(draws[, paste0("theta_", seq_len(K)), drop = FALSE])
  beta_matrix <- if (P > 0) {
    as.matrix(draws[, paste0("beta_", seq_len(P)), drop = FALSE])
  } else {
    NULL
  }
  sigma2 <- draws$sigma2
  sigma <- sqrt(sigma2)
  n_draws <- nrow(draws)

  time <- data$time
  event <- data$event
  group <- data$group
  X <- data$X
  n <- data$n

  log_time <- log(time)
  # linear predictor per observation per draw: eta[i, s] = theta[s, group[i]]
  eta <- t(theta_matrix[, group, drop = FALSE])
  if (P > 0) {
    eta <- eta + X %*% t(beta_matrix)
  }

  # pointwise log-likelihood matrix (n x n_draws)
  logpd <- matrix(NA_real_, nrow = n, ncol = n_draws)
  for (i in seq_len(n)) {
    llt <- (log_time[i] - eta[i, ]) / sigma
    if (event[i] == 1) {
      # log density of log-normal at (log t)
      logpd[i, ] <- stats::dnorm(llt, log = TRUE) - log(sigma) - log_time[i]
    } else {
      logpd[i, ] <- pnorm(llt, lower.tail = FALSE, log.p = TRUE)
    }
  }

  elpd_loo_onestep <- apply(logpd, 1, function(row) {
    # log mean exp
    m <- max(row)
    m + log(mean(exp(row - m)))
  })
  lpd <- sum(elpd_loo_onestep)

  # p_waic as sample variance of pointwise log-likelihood
  p_waic_vec <- apply(logpd, 1, function(row) stats::var(row) * (n_draws - 1) / n_draws)
  p_waic <- sum(p_waic_vec)

  waic <- -2 * (lpd - p_waic)
  # standard error across observations (per loo convention: se of sum of lpd)
  se_points <- sqrt(n * stats::var(elpd_loo_onestep))

  elpd_waic <- lpd - p_waic

  pointwise <- cbind(lpd = elpd_loo_onestep, p_waic = p_waic_vec, waic = waic)

  list(
    waic = waic,
    se_waic = se_points,
    lpd = lpd,
    p_waic = p_waic,
    elpd_waic = elpd_waic,
    pointwise = pointwise
  )
}

#' Compare multiple exnex_surv fits by WAIC
#'
#' Computes [compute_waic()] for each supplied fit and reports them in a
#' single `data.frame`, sorted by ascending WAIC.
#'
#' @param ... Two or more fitted `exnex_surv` objects.
#' @param digit Number of decimal places for the reported statistics.
#'
#' @return A `data.frame` with one row per model and columns `model`, `waic`,
#'   `se_waic`, `lpd`, `p_waic`, `elpd_waic`.
#' @export
compare_waic <- function(..., digit = 2) {
  fits <- list(...)
  if (length(fits) < 2) {
    stop("`compare_waic()` requires at least two fits.", call. = FALSE)
  }
  for (f in fits) checkmate::assert_class(f, "exnex_surv")

  fit_names <- names(fits)
  if (is.null(fit_names)) fit_names <- paste0("fit", seq_along(fits))

  out <- lapply(seq_along(fits), function(i) {
    w <- compute_waic(fits[[i]])
    data.frame(
      model = fit_names[i],
      waic = w$waic,
      se_waic = w$se_waic,
      lpd = w$lpd,
      p_waic = w$p_waic,
      elpd_waic = w$elpd_waic
    )
  })

  out <- do.call(rbind, out)
  out <- out[order(out$waic), ]
  rownames(out) <- NULL
  out[, c("waic", "se_waic")] <- round(out[, c("waic", "se_waic")], digit)
  out
}
