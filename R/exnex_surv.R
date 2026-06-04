#' Fit Bayesian EXNEX Survival Models
#'
#' Fits a Bayesian Exchangeable Non-Exchangeable (EXNEX) hierarchical model
#' for right-censored survival data in basket trials.
#'
#' The first predictor variable in the model formula is interpreted as the
#' basket/group assignment. Additional variables are included as covariates
#' in a linear predictor on the log-survival scale.
#'
#' @param x An object containing the predictors (subgroup assignment).
#'   Can be a data frame or a formula.
#' @param ... Additional arguments passed to methods.
#' @param formula A model formula with structure: `Surv(time, event) ~ group + covariates`.
#'   The first RHS variable is the group/basket assignment. Additional variables are covariates.
#'
#' @return An object of class `exnex_surv` containing posterior samples and metadata.
#'
#' @examples
#' # Fit model with group effects only
#' # fit <- exnex_surv(Surv(time, event) ~ group, data = trial_data)
#'
#' # Fit model with group + covariates
#' # fit <- exnex_surv(Surv(time, event) ~ group + age + baseline, data = trial_data)
#'
#' @export
exnex_surv <- function(x, ...) {
  UseMethod("exnex_surv")
}

#' @export
#' @rdname exnex_surv
exnex_surv.default <- function(x, ...) {
  stop(
    "`exnex_surv()` is not defined for a '",
    class(x)[1],
    "'.",
    call. = FALSE
  )
}

#' @param data A data frame containing the variables in the formula.
#' @param priors A list of hyperparameters. Currently a placeholder.
#' @param iter Total number of MCMC iterations. Default is 2000.
#' @param warmup Number of warmup iterations to discard. Default is 1000.
#'   Posterior samples will have (iter - warmup) rows.
#' @param chains Number of independent MCMC chains. Default is 1.
#' @param group_col Name of the column that represents the basket/group assignment.
#'   This variable will be treated as the group index, separate from covariates.
#'   If NULL (default), assumes the first RHS variable in the formula is the group.
#' @param seed Random seed for reproducibility (optional).
#'
#' @export
#' @rdname exnex_surv
exnex_surv.formula <- function(
  formula,
  data,
  priors = list(),
  iter = 2000,
  warmup = 1000,
  chains = 1,
  group_col = NULL,
  seed = NULL,
  ...
) {
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data, min.rows = 1, min.cols = 2)
  checkmate::assert_list(priors)
  checkmate::assert_int(iter, lower = 1)
  checkmate::assert_int(warmup, lower = 0)
  checkmate::assert_int(chains, lower = 1)
  checkmate::assert_character(group_col, len = 1, null.ok = TRUE)
  checkmate::assert_int(seed, lower = 1, upper = 2147483647, null.ok = TRUE)

  if (warmup >= iter) {
    stop(
      "`warmup` (",
      warmup,
      ") must be less than `iter` (",
      iter,
      ").",
      call. = FALSE
    )
  }

  if (!is.null(group_col) && !group_col %in% colnames(data)) {
    stop(
      "Column '",
      group_col,
      "' not found in data.",
      call. = FALSE
    )
  }

  if (is.null(group_col)) {
    rhs_vars <- all.vars(formula[[3]])
    if (length(rhs_vars) == 0) {
      stop(
        "Formula must have at least one RHS variable (the group).",
        call. = FALSE
      )
    }
    group_col <- rhs_vars[1]
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  processed <- hardhat::mold(formula, data)

  exnex_surv_bridge(
    processed = processed,
    priors = priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    group_col = group_col,
    original_data = data
  )
}

#' @param y A Surv object or matrix containing outcome (time and event status).
#' @param chains Number of independent MCMC chains. Default is 1.
#' @param group_col Name of the column in `x` that represents the basket/group assignment.
#'   If NULL (default), assumes the first column in x is the group.
#' @param seed Random seed for reproducibility (optional).
#'
#' @export
#' @rdname exnex_surv
exnex_surv.data.frame <- function(
  x,
  y,
  priors = list(),
  iter = 2000,
  warmup = 1000,
  chains = 1,
  group_col = NULL,
  seed = NULL,
  ...
) {
  checkmate::assert_data_frame(x, min.rows = 1, min.cols = 1)

  if (!inherits(y, "Surv") && !is.data.frame(y) && !is.matrix(y)) {
    stop(
      "`y` must be a Surv object, data frame, or matrix. ",
      "Got ",
      class(y)[1],
      ".",
      call. = FALSE
    )
  }

  checkmate::assert_list(priors)
  checkmate::assert_int(iter, lower = 1)
  checkmate::assert_int(warmup, lower = 0)
  checkmate::assert_int(chains, lower = 1)
  checkmate::assert_character(group_col, len = 1, null.ok = TRUE)
  checkmate::assert_int(seed, lower = 1, upper = 2147483647, null.ok = TRUE)

  n_y <- if (inherits(y, "Surv")) nrow(y) else nrow(y)

  if (nrow(x) != n_y) {
    stop(
      "Number of rows in predictors (",
      nrow(x),
      ") does not match outcomes (",
      n_y,
      ").",
      call. = FALSE
    )
  }

  if (warmup >= iter) {
    stop(
      "`warmup` (",
      warmup,
      ") must be less than `iter` (",
      iter,
      ").",
      call. = FALSE
    )
  }

  if (!is.null(group_col) && !group_col %in% colnames(x)) {
    stop(
      "Column '",
      group_col,
      "' not found in predictors.",
      call. = FALSE
    )
  }

  if (is.null(group_col)) {
    group_col <- colnames(x)[1]
  }

  if (!is.null(seed)) {
    set.seed(seed)
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

  exnex_surv_bridge(
    processed = processed,
    priors = priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    group_col = group_col,
    original_data = original_x
  )
}
