#' Fit Bayesian Pooling Survival Models
#'
#' Fits Bayesian models for right-censored log-normal survival data in basket
#' trials using a data-augmented Gibbs sampler. The `pooling` argument selects
#' the structural assumption on the basket effects:
#'
#' For patient \eqn{i}, the model is
#' \deqn{\log T_i = \theta_{g[i]} + X_i^{\mathsf T}\beta + \varepsilon_i,\quad
#' \varepsilon_i\sim\mathcal N(0,\sigma^2),}
#' and the argument \code{pooling} selects the prior structure of the basket
#' effects \eqn{\theta_j}:
#' \describe{
#'   \item{\code{"exnex"} (default)}{latent EXNEX hierarchy,
#'     \deqn{\theta_j\mid Z_j\sim Z_j\,\mathcal N(\mu,\tau^2)
#'     +(1-Z_j)\,\mathcal N(m_{0j},v_{0j}),\quad
#'     Z_j\sim\mathrm{Bern}(p_{\mathrm{exch},j}),}
#'     so that each basket either borrows strength from the exchangeable
#'     component \eqn{\mathcal N(\mu,\tau^2)} or follows its own
#'     non-exchangeable prior \eqn{\mathcal N(m_{0j},v_{0j})}.}
#'   \item{\code{"complete"}}{complete pooling: a single shared basket effect
#'     \eqn{\theta_j=\theta} for all baskets with prior
#'     \eqn{\theta\sim\mathcal N(m_{\mu},v_{\mu})}.}
#'   \item{\code{"none"}}{no pooling: independent basket effects,
#'     \eqn{\theta_j\sim\mathcal N(m_{0j},v_{0j})} for \eqn{j=1,\dots,K}.}
#' }
#' The first predictor variable in the model formula is interpreted as the
#' basket/group assignment; additional variables are covariates in the linear
#' predictor on the log-survival scale.
#'
#' Censored event times are handled by data augmentation: a censored log-time is
#' imputed from its truncated-Normal conditional distribution (inverse-CDF in log
#' space) before the remaining parameters are updated with their conjugate full
#' conditionals. The sampler keeps only the posterior draws of \eqn{\theta_j},
#' \eqn{\beta}, and \eqn{\sigma^2}.
#'
#' @param x An object containing the predictors (subgroup assignment).
#'   Can be a data frame or a formula.
#' @param ... Additional arguments passed to methods.
#' @param formula A model formula with structure: `Surv(time, event) ~ group + covariates`.
#'   The first RHS variable is the group/basket assignment. Additional variables are covariates.
#' @param pooling Pooling mode for the basket effects: one of
#'   \code{"exnex"} (default), \code{"complete"}, or \code{"none"} (see
#'   Description).
#'
#' @return An object of class `pooling_surv` containing posterior samples and metadata.
#'   Components include \code{pooling} (the fitted variant), \code{draws} (a data
#'   frame with columns \code{theta_1}, ..., \code{theta_K}, \code{beta_1}, ...,
#'   \code{beta_P}, and \code{sigma2}; under \code{pooling = "complete"} the
#'   \code{theta_*} columns hold identical draws from the shared effect),
#'   \code{data} (the processed data: \code{time}, \code{event},
#'   \code{group}, \code{X}, \code{n}, \code{n_groups}, \code{n_covariates},
#'   \code{cov_names}, \code{chain_seeds}), \code{priors} (the supplied priors),
#'   \code{resolved_priors} (defaults merged with overrides), \code{iter},
#'   \code{warmup}, \code{chains}, and \code{blueprint}. Use \code{summary()},
#'   \code{print()}, and \code{plot()} to inspect it.
#'
#' @examples
#' # Small simulated dataset: three baskets, one covariate, ~25% censoring
#' set.seed(1)
#' n <- 90
#' group <- factor(rep(1:3, each = 30))
#' x1 <- rnorm(n)
#' eta <- rep(c(1.1, 1.6, 2.0), each = 30) + 0.5 * x1
#' true_time <- exp(eta + rnorm(n, 0, 0.6))
#' cens_time <- runif(n, 2, 9)
#' d <- data.frame(
#'   time = pmin(true_time, cens_time),
#'   event = as.integer(true_time <= cens_time),
#'   group = group,
#'   x1 = x1
#' )
#'
#' # Fit the EXNEX hierarchy with group effects and one covariate
#' fit <- pooling_surv(
#'   survival::Surv(time, event) ~ group + x1,
#'   data = d,
#'   priors = list(p_mix = 0.7),
#'   iter = 300, warmup = 150, chains = 1
#' )
#' print(fit, show_trace = FALSE)
#' summary(fit)
#'
#' # Complete and no pooling use the same data augmentation machinery
#' fit_complete <- pooling_surv(
#'   survival::Surv(time, event) ~ group + x1,
#'   data = d, pooling = "complete",
#'   iter = 300, warmup = 150, chains = 1
#' )
#' fit_none <- pooling_surv(
#'   survival::Surv(time, event) ~ group + x1,
#'   data = d, pooling = "none",
#'   iter = 300, warmup = 150, chains = 1
#' )
#'
#' @export
pooling_surv <- function(x, ...) {
  UseMethod("pooling_surv")
}

#' @export
#' @rdname pooling_surv
pooling_surv.default <- function(x, ...) {
  stop(
    "`pooling_surv()` is not defined for a '",
    class(x)[1],
    "'.",
    call. = FALSE
  )
}

#' Shared validation for common pooling_surv arguments
#' @keywords internal
validate_args <- function(priors, iter, warmup, chains, parallel_chains, group_col, seed) {
  checkmate::assert_list(priors)
  checkmate::assert_int(iter, lower = 1)
  checkmate::assert_int(warmup, lower = 0)
  checkmate::assert_int(chains, lower = 1)
  checkmate::assert_int(parallel_chains, lower = 1)
  if (parallel_chains > chains) {
    stop("`parallel_chains` (", parallel_chains, ") must be less than or equal to `chains` (", chains, ").", call. = FALSE)
  }
  checkmate::assert_character(group_col, len = 1, null.ok = TRUE)
  checkmate::assert_int(seed, lower = 1, upper = 2147483647, null.ok = TRUE)
  if (warmup >= iter) {
    stop("`warmup` (", warmup, ") must be less than `iter` (", iter, ").", call. = FALSE)
  }
}

#' Validate the pooling argument
#' @keywords internal
validate_pooling <- function(pooling) {
  checkmate::assert_choice(
    pooling,
    c("exnex", "complete", "none"),
    .var.name = "pooling"
  )
}

#' @param data A data frame containing the variables in the formula.
#' @param priors Optional named list of prior hyperparameters. Supported fields:
#'   \code{a_sigma}, \code{b_sigma} (inverse-Gamma shape and scale for the
#'   residual variance \eqn{\sigma^2}); \code{a_tau}, \code{b_tau}
#'   (inverse-Gamma for the between-basket variance \eqn{\tau^2}.
#'   Only used with \code{pooling = "exnex"}); \code{p_mix} (EXNEX mixture
#'   weight, strictly between 0 and 1. Only used with
#'   \code{pooling = "exnex"}); \code{m_mu}, \code{v_mu} (prior mean and
#'   variance of the exchangeable center \eqn{\mu}; with
#'   \code{pooling = "complete"} they act as the prior of the shared basket
#'   effect); \code{m_nex}, \code{v_nex} (prior mean and variance of the
#'   nonexchangeable component; with \code{pooling = "none"} they are the
#'   priors of the independent basket effects); \code{v_beta}
#'   (variance of the regression-coefficient prior). \code{p_mix}, \code{m_nex},
#'   and \code{v_nex} each accept either a scalar, replicated across baskets, or
#'   a numeric vector of length K with one value per basket, matching the
#'   basket-specific notation of the model. All other fields are scalars. Absent
#'   fields keep the defaults
#'   (inverse-Gamma(2,2) for the variances, \code{p_mix = 0.5},
#'   \code{m_mu = 0}, \code{v_mu = 1e4}, \code{m_nex = 0}, \code{v_nex = 1e4},
#'   \code{v_beta = 1e4}); unknown fields are ignored. Priors that are
#'   irrelevant for the selected \code{pooling} mode are ignored with a warning.
#' @param iter Total number of MCMC iterations. Default is 2000.
#' @param warmup Number of warmup iterations to discard. Default is 1000.
#'   Posterior samples will have (iter - warmup) rows.
#' @param chains Number of independent MCMC chains. Default is 1.
#' @param parallel_chains Number of chains to run in parallel at the R level.
#'   Must be between 1 and `chains`. Default is 1 (sequential chain execution).
#' @param group_col Name of the column that represents the basket/group assignment.
#'   This variable will be treated as the group index, separate from covariates.
#'   If NULL (default), assumes the first RHS variable in the formula is the group.
#' @param seed Random seed for reproducibility (optional).
#'
#' @export
#' @rdname pooling_surv
pooling_surv.formula <- function(
  formula,
  data,
  pooling = "exnex",
  priors = list(),
  iter = 2000,
  warmup = 1000,
  chains = 1,
  parallel_chains = 1,
  group_col = NULL,
  seed = NULL,
  ...
) {
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data, min.rows = 1, min.cols = 2)
  validate_pooling(pooling)
  validate_args(priors, iter, warmup, chains, parallel_chains, group_col, seed)

  if (!is.null(group_col) && !group_col %in% colnames(data)) {
    stop("Column '", group_col, "' not found in data.", call. = FALSE)
  }

  if (is.null(group_col)) {
    rhs_vars <- all.vars(formula[[3]])
    if (length(rhs_vars) == 0) {
      stop("Formula must have at least one RHS variable (the group).", call. = FALSE)
    }
    group_col <- rhs_vars[1]
    if (!group_col %in% colnames(data)) {
      stop("Column '", group_col, "' not found in data.", call. = FALSE)
    }
  }

  processed <- hardhat::mold(formula, data)

  pooling_surv_bridge(
    processed = processed,
    pooling = pooling,
    priors = priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    parallel_chains = parallel_chains,
    group_col = group_col,
    original_data = data,
    seed = seed
  )
}

#' @param y A Surv object or matrix containing outcome (time and event status).
#' @param chains Number of independent MCMC chains. Default is 1.
#' @param parallel_chains Number of chains to run in parallel at the R level.
#'   Must be between 1 and `chains`. Default is 1 (sequential chain execution).
#' @param group_col Name of the column in `x` that represents the basket/group assignment.
#'   If NULL (default), assumes the first column in x is the group.
#' @param seed Random seed for reproducibility (optional).
#'
#' @export
#' @rdname pooling_surv
pooling_surv.data.frame <- function(
  x,
  y,
  pooling = "exnex",
  priors = list(),
  iter = 2000,
  warmup = 1000,
  chains = 1,
  parallel_chains = 1,
  group_col = NULL,
  seed = NULL,
  ...
) {
  checkmate::assert_data_frame(x, min.rows = 1, min.cols = 1)

  if (!inherits(y, "Surv") && !is.data.frame(y) && !is.matrix(y)) {
    stop("`y` must be a Surv object, data frame, or matrix. Got ", class(y)[1], ".", call. = FALSE)
  }

  validate_pooling(pooling)
  validate_args(priors, iter, warmup, chains, parallel_chains, group_col, seed)

  n_y <- nrow(y)

  if (nrow(x) != n_y) {
    stop("Number of rows in predictors (", nrow(x), ") does not match outcomes (", n_y, ").", call. = FALSE)
  }

  if (!is.null(group_col) && !group_col %in% colnames(x)) {
    stop("Column '", group_col, "' not found in predictors.", call. = FALSE)
  }

  if (is.null(group_col)) {
    group_col <- colnames(x)[1]
  }

  original_x <- x

  if (inherits(y, "Surv") || is.matrix(y)) {
    y_df <- data.frame(
      time = as.numeric(y[, 1]),
      event = as.numeric(y[, 2])
    )
  } else {
    y_df <- y
  }

  processed <- hardhat::mold(x, y_df)

  pooling_surv_bridge(
    processed = processed,
    pooling = pooling,
    priors = priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    parallel_chains = parallel_chains,
    group_col = group_col,
    original_data = original_x,
    seed = seed
  )
}
