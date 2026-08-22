#' Summarize posterior draws from an exnex_surv fit
#'
#' @param object A fitted `exnex_surv` object.
#' @param probs Quantiles to report.
#' @param ... Unused.
#'
#' @return A data frame with posterior summaries.
#' @export
summary.exnex_surv <- function(object, probs = c(0.05, 0.5, 0.95), ...) {
  checkmate::assert_class(object, "exnex_surv")
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

  out
}

.as_mcmc_array_exnex_surv <- function(draws, chains) {
  checkmate::assert_data_frame(draws, min.rows = 1, min.cols = 1)
  checkmate::assert_int(chains, lower = 1)

  n_draws <- nrow(draws)
  if (n_draws %% chains != 0) {
    stop(
      "Number of draw rows must be divisible by number of chains.",
      call. = FALSE
    )
  }

  n_iter <- n_draws / chains
  arr <- array(
    as.matrix(draws),
    dim = c(as.integer(n_iter), chains, ncol(draws))
  )

  dimnames(arr) <- list(
    iteration = as.character(seq_len(as.integer(n_iter))),
    chain = paste0("chain_", seq_len(chains)),
    parameter = colnames(draws)
  )

  arr
}

#' Plot parameter traces from an exnex_surv fit
#'
#' Generates one bayesplot traceplot per parameter.
#'
#' @param x A fitted `exnex_surv` object.
#' @param parameters Optional character vector of parameter names to plot.
#' @param ask Should R pause between plots? Defaults to `interactive()`.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
plot.exnex_surv <- function(x, parameters = NULL, ask = interactive(), ...) {
  checkmate::assert_class(x, "exnex_surv")

  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop(
      "Package 'bayesplot' is required for plot.exnex_surv().",
      call. = FALSE
    )
  }

  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop(
      "Package 'viridisLite' is required for plot.exnex_surv().",
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

  arr <- .as_mcmc_array_exnex_surv(draws, chains = x$chains)

  for (param in parameters) {
    trace_plot <- bayesplot::mcmc_trace(arr, pars = param)
    print(trace_plot)
  }

  invisible(x)
}

#' Print an exnex_surv fit
#'
#' Prints a compact summary and optionally traceplots of posterior draws.
#'
#' @param x A fitted `exnex_surv` object.
#' @param show_trace Should traceplots be shown? Default is `TRUE`.
#' @param parameters Optional character vector of parameter names for traceplots.
#' @param max_parameters Maximum number of parameters to plot.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.exnex_surv <- function(
  x,
  show_trace = TRUE,
  parameters = NULL,
  max_parameters = Inf,
  ...
) {
  checkmate::assert_class(x, "exnex_surv")
  checkmate::assert_flag(show_trace)
  checkmate::assert_numeric(
    max_parameters,
    len = 1,
    lower = 1,
    any.missing = FALSE
  )

  cat("<exnex_surv model>\n")
  cat("Draws:", nrow(x$draws), "total post-warmup samples\n")
  cat("      ", nrow(x$draws) / x$chains, "post-warmup samples per chain\n")
  cat("Groups:", x$data$n_groups, "| Covariates:", x$data$n_covariates, "\n")
  cat(
    "MCMC:",
    "iter =",
    x$iter,
    ", warmup =",
    x$warmup,
    ", chains =",
    x$chains,
    "\n\n"
  )

  summ <- summary(x)
  print(
    summ[, c("parameter", "mean", "sd", "q05", "q50", "q95")],
    row.names = FALSE
  )

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
