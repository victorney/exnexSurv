#' Constructor for exnex_surv Objects
#'
#' Creates a formal S3 object containing posterior samples and metadata
#' from EXNEX model fitting.
#'
#' @param draws A data frame containing posterior samples.
#'   Columns are: theta_1, ..., theta_K, beta_1, ..., beta_P (if P > 0), sigma2.
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
#' @param iter Total MCMC iterations performed.
#' @param warmup Number of warmup iterations discarded.
#' @param chains Number of chains run.
#' @param blueprint The hardhat blueprint for the original formula/data structure.
#'
#' @return An `exnex_surv` object (S3 class).
#' @keywords internal
new_exnex_surv <- function(
  draws,
  data,
  priors,
  iter,
  warmup,
  chains,
  blueprint
) {
  checkmate::assert_data_frame(draws, min.rows = 1, min.cols = 1)
  checkmate::assert_list(data, min.len = 1)
  checkmate::assert_list(priors)
  checkmate::assert_int(iter, lower = 1)
  checkmate::assert_int(warmup, lower = 0)
  checkmate::assert_int(chains, lower = 1)

  checkmate::assert_numeric(data$time, any.missing = FALSE)
  checkmate::assert_numeric(
    data$event,
    len = length(data$time),
    any.missing = FALSE
  )
  checkmate::assert_numeric(
    data$group,
    len = length(data$time),
    any.missing = FALSE
  )
  checkmate::assert_int(data$n, lower = 1)
  checkmate::assert_int(data$n_groups, lower = 1)
  checkmate::assert_int(data$n_covariates, lower = 0)

  expected_rows <- iter - warmup
  if (nrow(draws) != expected_rows) {
    stop(
      "Number of draw rows (",
      nrow(draws),
      ") does not match ",
      "iter - warmup (",
      expected_rows,
      ").",
      call. = FALSE
    )
  }

  expected_cols <- data$n_groups + data$n_covariates + 1 # K + P + sigma2
  if (ncol(draws) != expected_cols) {
    stop(
      "Number of columns in draws (",
      ncol(draws),
      ") does not match ",
      "expected K + P + 1 = ",
      expected_cols,
      ".",
      call. = FALSE
    )
  }

  object <- hardhat::new_model(
    draws = draws,
    data = data,
    priors = priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    blueprint = blueprint,
    class = "exnex_surv"
  )

  return(object)
}
