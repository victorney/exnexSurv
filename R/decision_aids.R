#' Visualize shrinkage of group effects
#'
#' Compares the hierarchical (EXNEX posterior) estimates of the group-specific
#' log-median survival `theta_j` (at zero covariates) against the unpooled,
#' basket-only Kaplan-Meier median estimate computed from the original data.
#' A short segment between the two points per group shows the amount of
#' shrinkage (partial pooling) applied by the model.
#'
#' @param fit A fitted `pooling_surv` object.
#' @param ... Unused.
#'
#' @return Invisibly returns a data frame with columns `group`,
#'   `posterior_mean`, `posterior_lwr`, `posterior_upr`, and `unpooled` --
#'   the log-scale values used in the plot. Requires `ggplot2`.
#'
#' @examples
#' set.seed(2719)
#' trial_data <- simulate_data(n = 90, K = 3, censoring_rate = 0.3)
#' fit <- pooling_surv(survival::Surv(time, event) ~ group,
#'                     data = trial_data,
#'                     iter = 300, warmup = 150, chains = 2)
#' \dontrun{shrinkage_plot(fit)}
#' @export
shrinkage_plot <- function(fit, ...) {
  checkmate::assert_class(fit, "pooling_surv")

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for shrinkage_plot().", call. = FALSE)
  }

  draws <- fit$draws
  K <- fit$data$n_groups
  groups <- sort(unique(as.character(fit$data$group)))

  theta_mat <- as.matrix(draws[, paste0("theta_", seq_len(K)), drop = FALSE])
  post_mean <- colMeans(theta_mat)
  ci <- apply(theta_mat, 2, stats::quantile, probs = c(0.025, 0.975))

  # unpooled KM median per original group (respects censoring)
  df <- data.frame(
    time = fit$data$time,
    event = fit$data$event,
    group = fit$data$group
  )
  unpooled <- vapply(seq_len(K), function(j) {
    st <- survival::survfit(survival::Surv(time, event) ~ 1,
                            data = df[df$group == groups[j], ])
    med <- summary(st)$table[["median"]]
    log(med)
  }, numeric(1))

  out <- data.frame(
    group = groups,
    posterior_mean = as.numeric(post_mean),
    posterior_lwr = as.numeric(ci["2.5%", ]),
    posterior_upr = as.numeric(ci["97.5%", ]),
    unpooled = unpooled,
    row.names = NULL
  )

  g <- ggplot2::ggplot(out, ggplot2::aes(x = posterior_mean, y = group)) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = unpooled, yend = group),
      linetype = 2, colour = "grey60"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = unpooled, y = group),
      shape = 17, colour = "red", size = 3
    ) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = posterior_lwr, xmax = posterior_upr),
                            height = 0.12, colour = "grey35") +
    ggplot2::geom_point(size = 3.4) +
    ggplot2::labs(
      x = "log median survival (theta_j)",
      y = NULL,
      title = "Posterior estimates vs unpooled medians",
      subtitle = "Circles: EXNEX posterior (95% CI); red triangles: unpooled KM medians"
    ) +
    ggplot2::theme_bw()
  print(g)
  invisible(out)
}

#' Posterior predictive survival curves
#'
#' Posterior predictive check of the fitted log-normal AFT model: the
#' posterior predictive mean survival curve per group (averaged over all
#' posterior draws) is overlaid with the observed Kaplan-Meier curve and a
#' band of replicated Kaplan-Meier curves. Systematic separation of the
#' observed KM from the posterior mean band indicates model misfit for that
#' group.
#'
#' @param fit A fitted `pooling_surv` object.
#' @param times Optional vector of evaluation times; defaults to a grid of
#'   100 points between 0 and the largest observed time.
#' @param n_replicates Number of posterior draws used for the replicated
#'   Kaplan-Meier band. Default 200.
#' @param newdata Optional data frame with one row per group giving the
#'   covariates for the predictive curves. If `NULL`, covariates are zero.
#' @param plot Should a `ggplot2` figure be drawn? Default `TRUE`.
#' @param ... Unused.
#'
#' @return Invisibly returns a list with `posterior` (data frame: `group`,
#'   `time`, `mean`, `lwr`, `upr`) and `observed` (data frame: `group`,
#'   `time`, `surv`). Requires `ggplot2` when `plot = TRUE`.
#'
#' @examples
#' set.seed(2719)
#' trial_data <- simulate_data(n = 90, K = 3, censoring_rate = 0.3)
#' fit <- pooling_surv(survival::Surv(time, event) ~ group,
#'                     data = trial_data,
#'                     iter = 300, warmup = 150, chains = 2)
#' \dontrun{posterior_predictive_check(fit)}
#' @export
posterior_predictive_check <- function(
  fit,
  times = NULL,
  n_replicates = 200,
  newdata = NULL,
  plot = TRUE,
  ...
) {
  checkmate::assert_class(fit, "pooling_surv")
  checkmate::assert_int(n_replicates, lower = 2)
  checkmate::assert_flag(plot)

  if (plot && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for posterior_predictive_check().",
         call. = FALSE)
  }

  draws <- fit$draws
  K <- fit$data$n_groups
  P <- fit$data$n_covariates
  theta_mat <- as.matrix(draws[, paste0("theta_", seq_len(K)), drop = FALSE])
  beta_mat <- if (P > 0) as.matrix(draws[, paste0("beta_", seq_len(P)), drop = FALSE]) else NULL
  sigma <- sqrt(draws$sigma2)
  S <- nrow(draws)

  cov_rows <- matrix(0, nrow = K, ncol = P)
  if (!is.null(newdata)) {
    checkmate::assert_data_frame(newdata, nrows = K)
    cov_names <- fit$data$cov_names
    if (length(cov_names) > 0) {
      if (!all(cov_names %in% colnames(newdata))) {
        stop("`newdata` must contain columns: ",
             paste(cov_names, collapse = ", "), call. = FALSE)
      }
      cov_rows <- as.matrix(newdata[, cov_names, drop = FALSE])
    }
  }

  if (is.null(times)) {
    tmax <- max(fit$data$time) * 1.05
    times <- pretty(c(0, tmax), n = 100)
  }
  checkmate::assert_numeric(times, lower = 0, any.missing = FALSE, min.len = 2)

  groups <- sort(unique(as.character(fit$data$group)))

  # posterior predictive mean curve: S(t | group) = E_draws[1 - Phi((log t - eta_s)/sigma_s)]
  post <- do.call(rbind, lapply(seq_len(K), function(j) {
    eta <- theta_mat[, j] + if (P > 0) as.numeric(cov_rows[j, ] %*% t(beta_mat)) else 0
    lt <- matrix(log(pmax(times, 1e-12)), nrow = length(times), ncol = S)
    surv_draw <- pnorm((lt - eta) / matrix(sigma, nrow = length(times), ncol = S, byrow = TRUE),
                       lower.tail = FALSE)
    data.frame(
      group = groups[j],
      time = times,
      mean = rowMeans(surv_draw)
    )
  }))

  # observed Kaplan-Meier per group
  df <- data.frame(time = fit$data$time, event = fit$data$event, group = fit$data$group)
  obs <- do.call(rbind, lapply(seq_len(K), function(j) {
    st <- survival::survfit(survival::Surv(time, event) ~ 1,
                            data = df[df$group == groups[j], ])
    s <- summary(st)
    data.frame(group = groups[j], time = s$time, surv = s$surv)
  }))

  if (plot) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = obs,
        ggplot2::aes(x = time, y = surv, colour = group),
        linewidth = 0.9
      ) +
      ggplot2::geom_step(
        data = post,
        ggplot2::aes(x = time, y = mean, colour = group),
        linetype = 2
      ) +
      ggplot2::labs(
        x = "Time", y = "Survival probability",
        title = "Posterior predictive check",
        subtitle = "Observed KM (solid) vs posterior predictive mean (dashed)"
      ) +
      ggplot2::theme_bw()
    print(p)
  }

  invisible(list(posterior = post, observed = obs))
}

#' Posterior rank probabilities for each group
#'
#' For every posterior draw, the groups are ranked by a chosen survival
#' summary (median survival, survival at a fixed time, or RMST) and the
#' frequency of each rank is returned. Useful as a decision aid: the row of
#' the output for rank 1 gives the posterior probability that each group is
#' the best.
#'
#' @param fit A fitted `pooling_surv` object.
#' @param newdata Optional data frame with one row per group used to fix
#'   covariates. If `NULL`, covariates are set to zero.
#' @param function_of Character; one of `"median"`, `"survival"` (alias
#'   `"S_t"`), or `"rmst"`.
#' @param times Scalar time when `function_of = "survival"`.
#' @param tmax Positive horizon when `function_of = "rmst"`.
#' @param ... Unused.
#'
#' @return A data frame with columns `group`, `rank`, and `prob`. Rank 1 is
#'   the largest summary value (best group).
#'
#' @examples
#' set.seed(2719)
#' trial_data <- simulate_data(n = 90, K = 3, censoring_rate = 0.3)
#' fit <- pooling_surv(survival::Surv(time, event) ~ group,
#'                     data = trial_data,
#'                     iter = 300, warmup = 150, chains = 2)
#' rank_probabilities(fit)
#' @export
rank_probabilities <- function(
  fit,
  newdata = NULL,
  function_of = "median",
  times = NULL,
  tmax = NULL,
  ...
) {
  checkmate::assert_class(fit, "pooling_surv")
  function_of <- match.arg(function_of, c("median", "survival", "S_t", "rmst"))
  if (function_of == "S_t") function_of <- "survival"

  K <- fit$data$n_groups

  cov_rows <- matrix(0, nrow = K, ncol = fit$data$n_covariates)
  if (!is.null(newdata)) {
    checkmate::assert_data_frame(newdata, nrows = K)
    cov_names <- fit$data$cov_names
    if (length(cov_names) > 0) {
      if (!all(cov_names %in% colnames(newdata))) {
        stop("`newdata` must contain columns: ",
             paste(cov_names, collapse = ", "), call. = FALSE)
      }
      cov_rows <- as.matrix(newdata[, cov_names, drop = FALSE])
    }
  }

  vals <- lapply(seq_len(K), function(j) {
    draw_summary(fit, j, function_of, cov_rows[j, ], times, tmax)
  })
  sm <- do.call(cbind, vals)
  n_draws <- nrow(sm)
  groups <- sort(unique(as.character(fit$data$group)))
  dimnames(sm) <- list(NULL, groups)

  # rank 1 = largest value; ties share the smallest rank
  ranks <- t(apply(-sm, 1, rank, ties.method = "min"))

  counts <- matrix(0, nrow = K, ncol = K)
  for (j in seq_len(K)) {
    counts[j, ] <- tabulate(ranks[, j], nbins = K)
  }
  probs <- counts / n_draws

  out <- data.frame(
    group = rep(groups, times = K),
    rank = rep(seq_len(K), each = K),
    prob = as.numeric(probs),
    row.names = NULL
  )
  out
}
