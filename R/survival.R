#' Survival curves from an exnex_surv fit
#'
#' Computes posterior survival curves \eqn{S(t) = \Pr(T > t)} from a fitted
#' `exnex_surv` model for a log-normal AFT specification:
#' \deqn{\log T_i = \theta_{g[i]} + X_i^\top\beta + \varepsilon,\qquad
#'   \varepsilon\sim\mathcal N(0,\sigma^2).}
#' For a fixed time \eqn{t} and linear predictor \eqn{\eta}, the survival
#' probability of the log-normal distribution is
#' \deqn{S(t) = 1 - \Phi\!\left(\frac{\log t - \eta}{\sigma}\right),}
#' where \eqn{\sigma^2} is the residual variance. Posterior draws of
#' \eqn{(\theta, \beta, \sigma^2)} are propagated through this expression to
#' obtain a full posterior distribution of \eqn{S(t)} at each time point.
#'
#' If `newdata` is supplied, each row is evaluated at its own group and
#' covariate values; otherwise the covariates are fixed at zero and the
#' first group is used (with a warning if more than one group exists).
#'
#' @param fit A fitted `exnex_surv` object.
#' @param newdata Optional data frame with columns matching the covariates of
#'   the model. If it contains a `group` column, that is used for the group
#'   index; otherwise the first group is used for all rows.
#' @param times Optional numeric vector of times at which to evaluate the
#'   curves. If `NULL`, a sensible grid is built from the observed follow-up
#'   times.
#' @param level Credible-interval level (default `0.95`).
#' @param ... Unused.
#'
#' @return An object of class `survival_exnex` (also a `data.frame`) with long
#'   format columns: `time`, `median`, `lower`, `upper`, and `group`.
#' @export
survival_curves <- function(fit, newdata = NULL, times = NULL, level = 0.95, ...) {
  checkmate::assert_class(fit, "exnex_surv")
  checkmate::assert_data_frame(newdata, null.ok = TRUE)
  checkmate::assert_numeric(times, lower = 0, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_number(level, lower = 0, upper = 1, finite = TRUE)

  draws <- fit$draws
  data <- fit$data
  K <- data$n_groups
  P <- data$n_covariates

  theta_matrix <- as.matrix(draws[, paste0("theta_", seq_len(K)), drop = FALSE])
  beta_matrix <- if (P > 0) as.matrix(draws[, paste0("beta_", seq_len(P)), drop = FALSE]) else NULL
  sigma <- sqrt(draws$sigma2)
  n_draws <- nrow(draws)

  if (is.null(times)) {
    times <- seq(0, max(data$time), length.out = 100)
  }

  # Resolve evaluation grid for each newdata row
  if (is.null(newdata)) {
    nr <- 1L
    grp_i <- 1L
    cov_rows <- matrix(0, nrow = 1, ncol = P)
    labels <- "all"
  } else {
    nr <- nrow(newdata)
    grp_i <- rep(1L, nr)
    if ("group" %in% colnames(newdata)) {
      gn <- sort(unique(as.character(data$group)))
      grp_i <- match(as.character(newdata[["group"]]), gn)
      if (anyNA(grp_i)) stop("`newdata$group` contains unknown levels.", call. = FALSE)
    }
    cov_rows <- if (P > 0) {
      as.matrix(newdata[, data$cov_names, drop = FALSE])
    } else {
      matrix(numeric(0), nrow = nr, ncol = 0)
    }
    labels <- as.character(seq_len(nr))
  }

  # Build a long evaluation grid: for each newdata row, repeat all times
  n_t <- length(times)
  n_eval <- nr * n_t
  time_rep <- rep(times, times = nr)
  grp_rep <- rep(grp_i, each = n_t)
  lbl_rep <- rep(labels, each = n_t)
  cov_rep <- if (P > 0) {
    cov_rows[rep(seq_len(nr), each = n_t), , drop = FALSE]
  } else {
    matrix(numeric(0), n_eval, 0)
  }

  # Compute S(t) per evaluation point across all posterior draws
  S_matrix <- matrix(NA_real_, nrow = n_draws, ncol = n_eval)
  for (i in seq_len(n_eval)) {
    eta_draw <- theta_matrix[, grp_rep[i]]
    if (P > 0) {
      eta_draw <- eta_draw + as.numeric(cov_rep[i, ] %*% t(beta_matrix))
    }
    S_matrix[, i] <- pnorm(
      (log(time_rep[i]) - eta_draw) / sigma,
      lower.tail = FALSE
    )
  }

  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  med <- apply(S_matrix, 2, stats::quantile, probs = probs[2], na.rm = FALSE)
  lo  <- apply(S_matrix, 2, stats::quantile, probs = probs[1], na.rm = FALSE)
  hi  <- apply(S_matrix, 2, stats::quantile, probs = probs[3], na.rm = FALSE)

  out <- data.frame(
    time = time_rep,
    median = med,
    lower = lo,
    upper = hi,
    group = lbl_rep,
    stringsAsFactors = FALSE
  )

  attr(out, "level") <- level
  class(out) <- c("survival_exnex", "data.frame")
  out
}

#' Plot survival curves from an exnex_surv fit
#'
#' Draws the posterior median and credible band of the survival function for
#' each group. Requires `ggplot2`.
#'
#' @param x An object of class `survival_exnex` returned by
#'   [survival_curves()].
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
plot.survival_exnex <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.survival_exnex().",
         call. = FALSE)
  }

  g <- ggplot2::ggplot(x, ggplot2::aes(
    x = time, y = median,
    ymin = lower, ymax = upper,
    group = group
  )) +
    ggplot2::geom_ribbon(ggplot2::aes(fill = group), alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(colour = group)) +
    ggplot2::labs(
      x = "Time",
      y = "Survival probability S(t)",
      colour = "Group",
      fill  = "Group"
    ) +
    ggplot2::theme_bw()

  print(g)
  invisible(x)
}

#' Posterior median survival time
#'
#' For the log-normal AFT model the median survival time for a linear
#' predictor \eqn{\eta} is \eqn{t_{med} = \exp(\eta)}, since \eqn{S(t)=0.5}
#' when \eqn{\log t = \eta}. Posterior draws of \eqn{\theta} and \eqn{\beta}
#' therefore induce a posterior distribution of \eqn{t_{med}} whose quantiles
#' are reported.
#'
#' @param fit A fitted `exnex_surv` object.
#' @param newdata Optional data frame (one row gives one median).
#' @param level Credible-interval level (default `0.95`).
#' @param ... Unused.
#'
#' @return A `data.frame` with columns `group`, `median`, `lower`, `upper`.
#' @export
median_survival <- function(fit, newdata = NULL, level = 0.95, ...) {
  checkmate::assert_class(fit, "exnex_surv")
  checkmate::assert_data_frame(newdata, null.ok = TRUE)
  checkmate::assert_number(level, lower = 0, upper = 1, finite = TRUE)

  draws <- fit$draws
  data <- fit$data
  K <- data$n_groups
  P <- data$n_covariates
  theta_matrix <- as.matrix(draws[, paste0("theta_", seq_len(K)), drop = FALSE])
  beta_matrix <- if (P > 0) as.matrix(draws[, paste0("beta_", seq_len(P)), drop = FALSE]) else NULL
  n_draws <- nrow(draws)

  if (is.null(newdata)) {
    nr <- 1L; grp_i <- 1L; cov_rows <- matrix(0, 1, P); labels <- "all"
  } else {
    nr <- nrow(newdata)
    grp_i <- rep(1L, nr)
    if ("group" %in% colnames(newdata)) {
      gn <- sort(unique(as.character(fit$data$group)))
      grp_i <- match(as.character(newdata[["group"]]), gn)
      if (anyNA(grp_i)) stop("`newdata$group` contains unknown levels.", call. = FALSE)
    }
    cov_rows <- if (P > 0) {
      as.matrix(newdata[, data$cov_names, drop = FALSE])
    } else {
      matrix(numeric(0), nrow = nr, ncol = 0)
    }
    labels <- as.character(seq_len(nr))
  }

  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  out <- lapply(seq_len(nr), function(i) {
    eta_draw <- theta_matrix[, grp_i[i]]
    if (P > 0) eta_draw <- eta_draw + as.numeric(cov_rows[i, ] %*% t(beta_matrix))
    t_med <- exp(eta_draw)
    qu <- stats::quantile(t_med, probs = probs)
    data.frame(group = labels[i], median = qu[2], lower = qu[1], upper = qu[3])
  })
  do.call(rbind, out)
}

#' Restricted mean survival time (RMST) from an exnex_surv fit
#'
#' Computes the restricted mean survival time up to a horizon `tmax`:
#' \deqn{RMST(tmax) = \int_0^{tmax} S(t)\,dt.}
#' The integral is evaluated numerically (trapezoidal rule over a fine grid)
#' for every posterior draw, then summarised with posterior quantiles.
#'
#' @param fit A fitted `exnex_surv` object.
#' @param tmax A positive scalar horizon up to which the RMST is computed.
#'   Defaults to `max(fit$data$time)`.
#' @param newdata Optional data frame (one row gives one RMST).
#' @param level Credible-interval level (default `0.95`).
#' @param grid_points Resolution of the numerical quadrature (default `400`).
#' @param ... Unused.
#'
#' @return A `data.frame` with columns `group`, `rmst`, `lower`, `upper`.
#' @export
rmst <- function(fit, tmax, newdata = NULL, level = 0.95, grid_points = 400, ...) {
  checkmate::assert_class(fit, "exnex_surv")
  checkmate::assert_number(tmax, lower = 0, finite = TRUE)
  if (tmax <= 0) stop("`tmax` must be positive.", call. = FALSE)
  checkmate::assert_data_frame(newdata, null.ok = TRUE)
  checkmate::assert_number(level, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_int(grid_points, lower = 100)

  g_times <- seq(0, tmax, length.out = grid_points)
  n_g <- length(g_times)

  draws <- fit$draws
  data <- fit$data
  K <- data$n_groups
  P <- data$n_covariates
  theta_matrix <- as.matrix(draws[, paste0("theta_", seq_len(K)), drop = FALSE])
  beta_matrix <- if (P > 0) as.matrix(draws[, paste0("beta_", seq_len(P)), drop = FALSE]) else NULL
  sigma <- sqrt(draws$sigma2)
  n_draws <- nrow(draws)

  if (is.null(newdata)) {
    nr <- 1L; grp_i <- 1L; cov_rows <- matrix(0, 1, P); labels <- "all"
  } else {
    nr <- nrow(newdata)
    grp_i <- rep(1L, nr)
    if ("group" %in% colnames(newdata)) {
      gn <- sort(unique(as.character(data$group)))
      grp_i <- match(as.character(newdata[["group"]]), gn)
      if (anyNA(grp_i)) stop("`newdata$group` contains unknown levels.", call. = FALSE)
    }
    cov_rows <- if (P > 0) {
      as.matrix(newdata[, data$cov_names, drop = FALSE])
    } else {
      matrix(numeric(0), nrow = nr, ncol = 0)
    }
    labels <- as.character(seq_len(nr))
  }

  # For one (group, cov_row), compute RMST draws
  rmst_draws <- function(grp, cov_row) {
    eta_draw <- theta_matrix[, grp]
    if (P > 0) eta_draw <- eta_draw + as.numeric(cov_row %*% t(beta_matrix))
    out <- numeric(n_draws)
    for (s in seq_len(n_draws)) {
      S <- pnorm((log(g_times) - eta_draw[s]) / sigma[s], lower.tail = FALSE)
      out[s] <- sum(0.5 * (S[-1] + S[-n_g]) * diff(g_times))
    }
    out
  }

  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  out <- lapply(seq_len(nr), function(i) {
    rd <- rmst_draws(grp_i[i], cov_rows[i, , drop = FALSE])
    qu <- stats::quantile(rd, probs = probs)
    data.frame(group = labels[i], rmst = qu[2], lower = qu[1], upper = qu[3])
  })
  do.call(rbind, out)
}
