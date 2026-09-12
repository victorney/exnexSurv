#' Summarize posterior draws from an pooling_surv fit
#'
#' @param object A fitted `pooling_surv` object.
#' @param probs Quantiles to report.
#' @param ... Unused.
#'
#' @return A data frame with posterior summaries.
#' @export
summary.pooling_surv <- function(object, probs = c(0.05, 0.5, 0.95), ...) {
  checkmate::assert_class(object, "pooling_surv")
  checkmate::assert_numeric(probs, lower = 0, upper = 1, any.missing = FALSE)

  draws <- object$draws

  out <- data.frame(
    parameter = colnames(draws),
    mean = colMeans(draws),
    sd = apply(draws, 2, stats::sd),
    row.names = NULL
  )

  for (p in probs) {
    out[[paste0(
      "q",
      formatC(100 * p, width = 2, flag = "0", format = "f", digits = 0)
    )]] <-
      apply(draws, 2, stats::quantile, probs = p)
  }

  # attach convergence diagnostics (rhat / ESS) when they are available
  diag <- tryCatch(.convergence_diagnostics_silent(object), error = function(e) NULL)
  if (!is.null(diag)) {
    out <- merge(out, diag, by = "parameter", sort = FALSE)
  }

  out
}

#' Plot parameter traces from an pooling_surv fit
#'
#' Generates one bayesplot traceplot per parameter.
#'
#' @param x A fitted `pooling_surv` object.
#' @param parameters Optional character vector of parameter names to plot.
#' @param ask Should R pause between plots? Defaults to `interactive()`.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
plot.pooling_surv <- function(x, parameters = NULL, ask = interactive(), ...) {
  checkmate::assert_class(x, "pooling_surv")

  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop(
      "Package 'bayesplot' is required for plot.pooling_surv().",
      call. = FALSE
    )
  }

  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop(
      "Package 'viridisLite' is required for plot.pooling_surv().",
      call. = FALSE
    )
  }

  draws <- x$draws
  all_parameters <- colnames(draws)

  if (is.null(parameters)) {
    parameters <- all_parameters
  }

  checkmate::assert_character(parameters, any.missing = FALSE, min.len = 1)

  missing_parameters <- setdiff(parameters, all_parameters)
  if (length(missing_parameters) > 0) {
    stop(
      "Unknown parameter(s): ",
      paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }

  viridis_scheme <- substr(viridisLite::viridis(6), 1, 7)
  bayesplot::color_scheme_set(viridis_scheme)

  old_ask <- graphics::par("ask")
  on.exit(graphics::par(ask = old_ask), add = TRUE)
  graphics::par(ask = ask)

  n_draws <- nrow(draws)
  n_iter <- n_draws / x$chains
  arr <- array(as.matrix(draws), dim = c(as.integer(n_iter), x$chains, ncol(draws)))
  dimnames(arr) <- list(iteration = as.character(seq_len(as.integer(n_iter))), chain = paste0("chain_", seq_len(x$chains)), parameter = colnames(draws))

  for (param in parameters) {
    trace_plot <- bayesplot::mcmc_trace(arr, pars = param)
    print(trace_plot)
  }

  invisible(x)
}

#' Print an pooling_surv fit
#'
#' Prints a compact summary and optionally traceplots of posterior draws.
#'
#' @param x A fitted `pooling_surv` object.
#' @param show_trace Should traceplots be shown? Default is `TRUE`.
#' @param parameters Optional character vector of parameter names for traceplots.
#' @param max_parameters Maximum number of parameters to plot.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.pooling_surv <- function(
  x,
  show_trace = TRUE,
  parameters = NULL,
  max_parameters = Inf,
  ...
) {
  checkmate::assert_class(x, "pooling_surv")
  checkmate::assert_flag(show_trace)
  checkmate::assert_numeric(
    max_parameters,
    len = 1,
    lower = 1,
    any.missing = FALSE
  )

  cat("<pooling_surv model>\n")
  cat("Pooling:", x$pooling, "\n")
  cat("Draws:", nrow(x$draws), "total post-warmup samples\n")
  cat("      ", nrow(x$draws) / x$chains, "post-warmup samples per chain\n")
  cat("Groups:", x$data$n_groups, "| Covariates:", x$data$n_covariates, "\n")
  cat(
    "MCMC:",
    "iter =",
    x$iter,
    ", warmup =",
    x$warmup,
    ", thin =",
    if (is.null(x$thin)) 1L else x$thin,
    ", chains =",
    x$chains,
    "\n\n"
  )

  summ <- summary(x)
  show_cols <- intersect(
    c("parameter", "mean", "sd", "q05", "q50", "q95", "rhat", "ess_bulk", "ess_tail"),
    colnames(summ)
  )
  print(summ[, show_cols], row.names = FALSE)

  if ("rhat" %in% colnames(summ) && any(is.finite(summ$rhat))) {
    worst_rhat <- max(summ$rhat, na.rm = TRUE)
    min_ess <- min(c(summ$ess_bulk, summ$ess_tail), na.rm = TRUE)
    cat(
      "\nConvergence: max R-hat =",
      formatC(worst_rhat, digits = 3), "| min ESS =",
      formatC(min_ess, digits = 1, format = "f")
    )
    if (worst_rhat > 1.01 || min_ess < 100) {
      cat(" (consider longer chains or more warmup)")
    }
    cat("\n")
  }

  if (show_trace) {
    if (!requireNamespace("bayesplot", quietly = TRUE)) {
      cat(
        "\nTraceplots skipped: package 'bayesplot' is not installed.\n"
      )
    } else {
      if (is.null(parameters)) {
        parameters <- colnames(x$draws)
      }

      n_to_plot <- if (is.infinite(max_parameters)) {
        length(parameters)
      } else {
        min(length(parameters), as.integer(max_parameters))
      }

      parameters <- parameters[seq_len(n_to_plot)]

      cat("\nTraceplots (one parameter per panel):\n")
      plot(x, parameters = parameters, ask = interactive())
    }
  }

  invisible(x)
}
