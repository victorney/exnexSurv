#' Convergence diagnostics (R-hat and effective sample size)
#'
#' Computes R-hat, bulk effective sample size (ESS), and tail ESS for every
#' parameter of a fitted `pooling_surv` object, using the `posterior`
#' package. All three follow the definitions in Vehtari, Gelman, Simpson,
#' Carpenter & Buerkner (2021, "Rank-normalization, folding, and
#' localization"): the R-hat is split and rank-normalized, bulk ESS is
#' computed on rank-normalized draws, and tail ESS quantifies sampling
#' adequacy in the 5%/95% tails. Each chain is split in half before the
#' diagnostic is computed (so `chains = 1` still gets a meaningful R-hat).
#'
#' Rules of thumb: R-hat close to 1.0 (say < 1.01) and ESS in the hundreds
#' give reasonably stable summaries. Low ESS with a good R-hat means the
#' chain is slow to move but honest; both R-hat above ~1.01 and very low
#' ESS (e.g. below 100) are worth investigating before trusting the
#' posterior.
#'
#' @param fit A fitted `pooling_surv` object (result of `pooling_surv()`).
#' @param ... Unused.
#'
#' @return A data frame with one row per parameter and columns `parameter`,
#'   `rhat`, `ess_bulk`, and `ess_tail`. All values are `NA` if the
#'   `posterior` package is not installed.
#'
#' @examples
#' set.seed(2719)
#' trial_data <- simulate_data(n = 90, K = 3, censoring_rate = 0.3)
#' fit <- pooling_surv(survival::Surv(time, event) ~ group,
#'                     data = trial_data,
#'                     iter = 300, warmup = 150, chains = 2)
#' convergence_diagnostics(fit)
#' @export
convergence_diagnostics <- function(fit, ...) {
  checkmate::assert_class(fit, "pooling_surv")

  params <- colnames(fit$draws)
  out <- data.frame(
    parameter = params,
    rhat = NA_real_,
    ess_bulk = NA_real_,
    ess_tail = NA_real_,
    row.names = NULL
  )

  if (!requireNamespace("posterior", quietly = TRUE)) {
    warning(
      "Package 'posterior' is not installed; returning NAs. ",
      "Install it with install.packages('posterior').",
      call. = FALSE
    )
    return(out)
  }

  # per-variable diagnostics on a (iteration, chain, variable) draws_array;
  # some posterior versions mis-handle multi-variable batch calls, so each
  # variable is subsetted individually
  arr <- .draws_array_pooling(fit)
  da <- posterior::as_draws_array(arr)
  for (v in seq_along(params)) {
    d_v <- posterior::subset_draws(da, variable = params[v])
    out$rhat[v] <- as.numeric(posterior::rhat(d_v))
    out$ess_bulk[v] <- as.numeric(posterior::ess_bulk(d_v))
    out$ess_tail[v] <- as.numeric(posterior::ess_tail(d_v))
  }
  out
}

#' Reshape combined draws into a (iteration, chain, parameter) array
#'
#' Draws rows are stored blockwise: chain 1 first, then chain 2, ...
#' @keywords internal
.draws_array_pooling <- function(fit) {
  m <- as.matrix(fit$draws)
  S <- nrow(m) / fit$chains
  if (S != floor(S)) {
    stop("draws rows are not an integer multiple of chains; draws frame is corrupt.",
         call. = FALSE)
  }
  a <- array(m, dim = c(as.integer(S), fit$chains, ncol(m)))
  dimnames(a) <- list(NULL, NULL, colnames(m))
  a
}

#' Convergence diagnostics without warnings (used internally by summary)
#' @keywords internal
.convergence_diagnostics_silent <- function(fit) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    return(NULL)
  }
  suppressWarnings(convergence_diagnostics(fit))
}
