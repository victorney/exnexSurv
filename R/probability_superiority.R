#' Posterior probability that one group beats another
#'
#' Estimates the posterior probability that a certain scalar summary of the
#' survival distribution for group `a` is greater than the same summary for
#' group `b`. Supported summaries are the median survival time, the survival
#' probability at a fixed time `S(t_0)`, and the restricted mean survival time
#' (RMST). The function works draw-by-draw on the posterior, so the probability
#' is computed on the joint posterior of the two groups.
#'
#' @param fit A fitted `exnex_surv` object.
#' @param a An integer index (or group label) of the first group.
#' @param b An integer index (or group label) of the second group.
#' @param newdata Optional data frame with one row per group, used to fix
#'   covariates. If `NULL`, covariates are set to zero.
#' @param function_of Character; one of `"median"`, `"survival"` (alias
#'   `"S_t"`), or `"rmst"`. Selects which summary is compared.
#' @param times A scalar time at which to evaluate the survival function when
#'   `function_of = "survival"`.
#' @param tmax A positive horizon when `function_of = "rmst"`.
#' @param ... Unused.
#'
#' @return A named list with elements `prob` (posterior probability that the
#'   summary of group `a` exceeds that of group `b`), `summary` (which summary
#'   was used), and `groups` (the two groups compared).
#' @export
probability_superiority <- function(
  fit,
  a,
  b,
  function_of = "median",
  newdata = NULL,
  times = NULL,
  tmax = NULL,
  ...
) {
  checkmate::assert_class(fit, "exnex_surv")
  function_of <- match.arg(function_of, c("median", "survival", "S_t", "rmst"))
  if (function_of == "S_t") function_of <- "survival"

  K <- fit$data$n_groups
  if (is.character(a)) a <- match(a, sort(unique(as.character(fit$data$group))))
  if (is.character(b)) b <- match(b, sort(unique(as.character(fit$data$group))))
  checkmate::assert_int(a, lower = 1, upper = K)
  checkmate::assert_int(b, lower = 1, upper = K)
  if (a == b) stop("`a` and `b` must be different groups.", call. = FALSE)

  if (function_of == "survival") {
    if (is.null(times)) stop("`times` is required when function_of = \"survival\".",
                             call. = FALSE)
    checkmate::assert_number(times, lower = 0, finite = TRUE)
  }
  if (function_of == "rmst") {
    if (is.null(tmax)) stop("`tmax` is required when function_of = \"rmst\".",
                            call. = FALSE)
    checkmate::assert_number(tmax, lower = 0, finite = TRUE)
  }

  # Resolve per-group covariate rows from `newdata`
  cov_row_a <- cov_row_b <- numeric(fit$data$n_covariates)
  if (!is.null(newdata)) {
    checkmate::assert_data_frame(newdata)
    if (nrow(newdata) != 2L) {
      stop("`newdata` must have exactly 2 rows: row 1 for group '",
           a, "', row 2 for group '", b, "'.", call. = FALSE)
    }
    # If a `group` column is present, use it to validate the expected groups.
    if ("group" %in% colnames(newdata)) {
      gn <- sort(unique(as.character(fit$data$group)))
      got <- as.character(newdata[["group"]])
      if (got[1] != gn[a] || got[2] != gn[b]) {
        stop("`newdata$group` must be c('", gn[a], "', '", gn[b],
             "'), got c('", got[1], "', '", got[2], "').", call. = FALSE)
      }
    }
    cov_names <- fit$data$cov_names
    if (length(cov_names) > 0) {
      if (!all(cov_names %in% colnames(newdata))) {
        stop("`newdata` must contain columns: ", paste(cov_names, collapse = ", "), call. = FALSE)
      }
      cov_row_a <- as.numeric(newdata[1, cov_names, drop = FALSE])
      cov_row_b <- as.numeric(newdata[2, cov_names, drop = FALSE])
    }
  }

  summary_a <- draw_summary(fit, a, function_of, cov_row_a, times, tmax)
  summary_b <- draw_summary(fit, b, function_of, cov_row_b, times, tmax)

  prob <- mean(summary_a > summary_b)

  list(
    prob = prob,
    summary = function_of,
    groups = c(a = a, b = b)
  )
}

# internal: posterior draws of a scalar survival summary for one group
draw_summary <- function(fit, group_idx, function_of, cov_row, times, tmax) {
  draws <- fit$draws
  K <- fit$data$n_groups
  P <- fit$data$n_covariates
  theta_matrix <- as.matrix(draws[, paste0("theta_", seq_len(K)), drop = FALSE])
  beta_matrix <- if (P > 0) {
    as.matrix(draws[, paste0("beta_", seq_len(P)), drop = FALSE])
  } else NULL
  sigma2 <- draws$sigma2
  n_draws <- nrow(draws)

  eta_draw <- theta_matrix[, group_idx]
  if (P > 0) eta_draw <- eta_draw + as.numeric(cov_row %*% t(beta_matrix))

  if (function_of == "median") {
    return(exp(eta_draw))
  }
  if (function_of == "survival") {
    return(pnorm((log(times) - eta_draw) / sqrt(sigma2), lower.tail = FALSE))
  }
  if (function_of == "rmst") {
    g_times <- seq(0, tmax, length.out = 400)
    out <- numeric(n_draws)
    for (s in seq_len(n_draws)) {
      S <- pnorm((log(g_times) - eta_draw[s]) / sqrt(sigma2[s]),
                 lower.tail = FALSE)
      out[s] <- sum(0.5 * (S[-1] + S[-length(S)]) * diff(g_times))
    }
    return(out)
  }
  stop("Unknown summary function.", call. = FALSE)
}

