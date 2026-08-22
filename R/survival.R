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
#' covariate values; otherwise the covariate values are fixed at zero and the
#' first group is used (with a warning if more than one group exists).
#'
#' @param fit A fitted `exnex_surv` object.
#' @param newdata Optional data frame with columns matching the covariates of
#'   the model (and optionally the group column). Each row yields one curve.
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

  theta_cols <- paste0("theta_", seq_len(K))
  beta_cols <- if (P > 0) paste0("beta_", seq_len(P)) else character(0)
  sigma2 <- draws$sigma2
  n_draws <- nrow(draws)

  theta_matrix <- as.matrix(draws[, theta_cols, drop = FALSE])
  beta_matrix <- if (P > 0) as.matrix(draws[, beta_cols, drop = FALSE]) else NULL

  # Resolve evaluation frame
  if (is.null(newdata)) {
    group_idx <- 1L
    cov_row <- numeric(P)
    label <- "all"
  } else {
    nr <- nrow(newdata)
    group_idx <- rep(1L, nr)
    if ("group" %in% colnames(newdata)) {
      g <- newdata[["group"]]
      if (is.factor(g)) g <- as.character(g)
      gn <- sort(unique(as.character(fit$data$group)))
      group_idx <- match(g, gn)
      if (anyNA(group_idx)) {
        stop("`newdata$group` contains levels not present in the fitted data.",
             call. = FALSE)
      }
    }
    if (P > 0) {
      cov_names <- data$cov_names
      if (!all(cov_names %in% colnames(newdata))) {
        stop("`newdata` must contain columns: ", paste(cov_names, collapse = ", "),
             call. = FALSE)
      }
      cov_row <- as.matrix(newdata[, cov_names, drop = FALSE])
    } else {
      cov_row <- matrix(numeric(0), nrow = nr, ncol = 0)
    }
    label <- as.character(seq_len(nr))
  }

  if (is.null(times)) {
    if (is.null(newdata)) {
      t0 <- 0
      tmax <- max(data$time)
    } else {
      t0 <- min(0, data$time)
      tmax <- max(data$time)
    }
    times <- seq(t0, tmax, length.out = 100)
  }

  # Broadcast evaluation points
  if (is.null(newdata)) {
    n_eval <- length(times)
    grp <- rep(1L, n_eval)
    cov_use <- matrix(0, nrow = n_eval, ncol = P)
    lbl <- rep(label, n_eval)
  } else {
    nr <- nrow(newdata)
    n_eval <- nr * length(times)
    grp <- rep(group_idx, each = length(times))
    cov_use <- if (P > 0) {
      cov_row[rep(seq_len(nr), each = length(times)), , drop = FALSE]
    } else {
      matrix(numeric(0), nrow = n_eval, ncol = 0)
    }
    lbl <- rep(label, each = length(times))
  }

  # Compute survival probability per observation across all posterior draws.
  S_matrix <- matrix(NA_real_, nrow = n_draws, ncol = n_eval)
  for (i in seq_len(n_eval)) {
    eta_draw <- theta_matrix[, grp[i]]
    if (P > 0) {
      eta_draw <- eta_draw + as.numeric(cov_use[i, ] %*% t(beta_matrix))
    }
    S_matrix[, i] <- pnorm(
      (rep(log(times)[i], n_draws) - eta_draw) / sqrt(sigma2),
      lower.tail = FALSE
    )
  }

  # Summarize S(t) per time point
  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  med <- apply(S_matrix, 2, stats::quantile, probs = probs[2])
  lo <- apply(S_matrix, 2, stats::quantile, probs = probs[1])
  hi <- apply(S_matrix, 2, stats::quantile, probs = probs[3])

  out <- data.frame(
    time = times,
    median = med,
    lower = lo,
    upper = hi,
    group = lbl,
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

  g <- ggplot2::ggplot(x, ggplot2::aes(x = .data$time, y = .data$median,
                                       ymin = .data$lower, ymax = .data$upper,
                                       group = .data$group)) +
    ggplot2::geom_ribbon(ggplot2::aes(fill = .data$group), alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(colour = .data$group)) +
    ggplot2::labs(
      x = "Time",
      y = "Survival probability S(t)",
      colour = "Group",
      fill = "Group"
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
  theta_cols <- paste0("theta_", seq_len(K))
  beta_cols <- if (P > 0) paste0("beta_", seq_len(P)) else character(0)

  theta_matrix <- as.matrix(draws[, theta_cols, drop = FALSE])
  beta_matrix <- if (P > 0) as.matrix(draws[, beta_cols, drop = FALSE]) else NULL
  n_draws <- nrow(draws)

  if (is.null(newdata)) {
    nr <- 1L
    group_idx <- 1L
    cov_row <- matrix(0, nrow = 1, ncol = P)
    label <- "all"
  } else {
    nr <- nrow(newdata)
    group_idx <- rep(1L, nr)
    if ("group" %in% colnames(newdata)) {
      gn <- sort(unique(as.character(data$group)))
      group_idx <- match(as.character(newdata[["group"]]), gn)
      if (anyNA(group_idx)) {
        stop("`newdata$group` contains levels not present in the fitted data.",
             call. = FALSE)
      }
    }
    cov_row <- if (P > 0) {
      as.matrix(newdata[, data$cov_names, drop = FALSE])
    } else {
      matrix(numeric(0), nrow = nr, ncol = 0)
    }
    label <- as.character(seq_len(nr))
  }

  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  out <- lapply(seq_len(nr), function(i) {
    eta_draw <- theta_matrix[, group_idx[i]]
    if (P > 0) {
      eta_draw <- eta_draw + as.numeric(cov_row[i, ] %*% t(beta_matrix))
    }
    t_med <- exp(eta_draw)
    qu <- stats::quantile(t_med, probs = probs)
    data.frame(group = label[i], median = qu[2], lower = qu[1], upper = qu[3])
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

  # Survival matrix for a single observation at this grid
  surv_at <- function(group_idx, cov_row) {
    draws <- fit$draws
    theta_matrix <- as.matrix(draws[, paste0("theta_", seq_len(fit$data$n_groups))])
    beta_matrix <- if (fit$data$n_covariates > 0) {
      as.matrix(draws[, paste0("beta_", seq_len(fit$data$n_covariates))])
    } else NULL
    sigma2 <- draws$sigma2
    n_draws <- nrow(draws)
    n_g <- length(g_times)
    S <- matrix(NA_real_, n_draws, n_g)
    eta_draw <- theta_matrix[, group_idx]
    if (!is.null(beta_matrix)) {
      eta_draw <- eta_draw + as.numeric(cov_row %*% t(beta_matrix))
    }
    for (j in seq_len(n_g)) {
      S[, j] <- pnorm(
        (rep(log(g_times[j]), n_draws) - eta_draw) / sqrt(sigma2),
        lower.tail = FALSE
      )
    }
    S
  }

  # Resolve evaluation rows
  if (is.null(newdata)) {
    nr <- 1L
    group_idx <- 1L
    cov_row <- matrix(0, 1, fit$data$n_covariates)
    label <- "all"
  } else {
    nr <- nrow(newdata)
    group_idx <- rep(1L, nr)
    if ("group" %in% colnames(newdata)) {
      gn <- sort(unique(as.character(fit$data$group)))
      group_idx <- match(as.character(newdata[["group"]]), gn)
      if (anyNA(group_idx)) {
        stop("`newdata$group` contains levels not present in the fitted data.",
             call. = FALSE)
      }
    }
    cov_row <- if (fit$data$n_covariates > 0) {
      as.matrix(newdata[, fit$data$cov_names, drop = FALSE])
    } else {
      matrix(0, nr, 0)
    }
    label <- as.character(seq_len(nr))
  }

  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  out <- lapply(seq_len(nr), function(i) {
    S <- surv_at(group_idx[i], cov_row[i, , drop = FALSE])
    grid_rmst <- apply(S, 1, function(srow) {
      sum(0.5 * (srow[-1] + srow[-ncol(S)]) * diff(g_times))
    })
    qu <- stats::quantile(grid_rmst, probs = probs)
    data.frame(group = label[i], rmst = qu[2], lower = qu[1], upper = qu[3])
  })

  do.call(rbind, out)
}
