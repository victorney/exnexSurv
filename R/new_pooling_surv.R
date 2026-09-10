#' Constructor for pooling_surv Objects
#'
#' Creates a formal S3 object containing posterior samples and metadata
#' from Bayesian pooling model fitting (EXNEX, complete pooling, or no
#' pooling).
#'
#' @param draws A data frame containing posterior samples.
#'   Columns are: theta_1, ..., theta_K, beta_1, ..., beta_P (if P > 0), sigma2.
#' @param pooling The fitted pooling variant: one of \code{"exnex"},
#'   \code{"complete"}, or \code{"none"}.
#' @param data A list containing:
#'   - time: Observed follow-up times
#'   - event: Event indicators (0/1)
#'   - group: Group assignments
#'   - X: Covariate matrix (can be empty)
#'   - n: Total number of observations
#'   - n_groups: Number of groups
#'   - n_covariates: Number of covariates
#'   - cov_names: Covariate column names
#' @param priors A list of prior specifications used for fitting.
#' @param resolved_priors A named list of the prior hyperparameters actually
#'   used (defaults merged with any overrides supplied in \code{priors}).
#' @param iter Total MCMC iterations performed.
#' @param warmup Number of warmup iterations discarded.
#' @param chains Number of chains run.
#' @param blueprint The hardhat blueprint for the original formula/data structure.
#'
#' @return An `pooling_surv` object (S3 class).
#' @keywords internal
new_pooling_surv <- function(
  draws,
  data,
  pooling,
  priors,
  resolved_priors = priors,
  iter,
  warmup,
  chains,
  blueprint
) {
  checkmate::assert_choice(pooling, c("exnex", "complete", "none"))
  checkmate::assert_data_frame(draws, min.rows = 1, min.cols = 1)
  checkmate::assert_list(data, min.len = 1)

  expected_rows <- (iter - warmup) * chains
  if (nrow(draws) != expected_rows) {
    stop("Number of draw rows (", nrow(draws), ") does not match iter - warmup (", expected_rows, ").", call. = FALSE)
  }

  expected_cols <- data$n_groups + data$n_covariates + 1
  if (ncol(draws) != expected_cols) {
    stop("Number of columns in draws (", ncol(draws), ") does not match expected K + P + 1 = ", expected_cols, ".", call. = FALSE)
  }

  object <- hardhat::new_model(
    draws = draws,
    data = data,
    pooling = pooling,
    priors = priors,
    resolved_priors = resolved_priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    blueprint = blueprint,
    class = "pooling_surv"
  )

  return(object)
}
