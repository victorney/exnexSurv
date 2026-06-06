#' Bridge connecting hardhat processed data to the C++ Gibbs Sampler
#'
#' Orchestrates data preparation, validation, and passes it to the C++ MCMC engine.
#'
#' The group variable (specified by group_col) is identified and reconstructed
#' from the hardhat-processed predictors. Additional variables are treated
#' as covariates in the linear predictor.
#'
#' @param processed A list produced by `hardhat::mold()`.
#' @param priors A list of hyperparameters (currently placeholder).
#' @param iter Total number of MCMC iterations.
#' @param warmup Number of warmup iterations.
#' @param chains Number of chains to run.
#' @param group_col Name of the original group column (before encoding).
#' @param original_data The original data frame (before hardhat processing).
#'
#' @return An `exnex_surv` object.
#' @keywords internal
exnex_surv_bridge <- function(
  processed,
  priors,
  iter,
  warmup,
  chains,
  group_col,
  original_data
) {
  predictors <- processed$predictors
  outcomes <- processed$outcomes

  # extract outcomes
  time_vec <- .extract_time_vector(outcomes)
  if (chains != 1) {
    stop(
      "Multi-chain support is not implemented yet. Use chains = 1.",
      call. = FALSE
    )
  }
  event_vec <- .extract_event_vector(outcomes)

  n <- length(time_vec)
  checkmate::assert_numeric(
    time_vec,
    len = n,
    finite = TRUE,
    any.missing = FALSE
  )
  checkmate::assert_numeric(event_vec, len = n, any.missing = FALSE)

  if (!all(time_vec > 0)) {
    stop("All survival times must be positive.", call. = FALSE)
  }

  if (!all(event_vec %in% c(0, 1))) {
    stop("Event status must be 0 (censored) or 1 (observed).", call. = FALSE)
  }

  # extract group assignment
  n_encoded_group_cols <- NA_integer_

  if (is.null(group_col)) {
    first_col <- predictors[[1]]

    if (is.numeric(first_col) && all(first_col %in% c(0, 1))) {
      tryCatch(
        {
          group_vec <- .extract_group_from_encoded(predictors)
          n_encoded_group_cols <- .count_group_columns(predictors)
        },
        error = function(e) {
          stop(
            "Failed to extract group from one-hot encoded predictors. ",
            "This may indicate hardhat encoding mismatch or invalid column names. ",
            "Original error: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
    } else if (is.factor(first_col)) {
      group_vec <- as.numeric(first_col)
      n_encoded_group_cols <- 1
    } else {
      stop(
        "First predictor must be either a factor or one-hot encoded numeric. ",
        "Got class: ",
        class(first_col)[1],
        call. = FALSE
      )
    }
  } else {
    if (!group_col %in% colnames(original_data)) {
      stop(
        "Column '",
        group_col,
        "' not found in original data.",
        call. = FALSE
      )
    }

    group_col_data <- original_data[[group_col]]

    if (is.list(group_col_data) && !inherits(group_col_data, "factor")) {
      stop(
        "group_col data is a list but not a factor. Got: ",
        class(group_col_data)[1],
        call. = FALSE
      )
    }

    group_vec <- as.numeric(as.factor(group_col_data))
    n_encoded_group_cols <- NA_integer_
  }

  if (any(is.na(group_vec))) {
    stop("Group variable contains missing values.", call. = FALSE)
  }

  n_groups <- max(group_vec)
  checkmate::assert_int(n_groups, lower = 1)

  # extract covariates
  if (is.null(group_col)) {
    if (ncol(predictors) > n_encoded_group_cols) {
      cols_to_remove <- 1:n_encoded_group_cols
      X_mat <- as.matrix(predictors[, -cols_to_remove, drop = FALSE])
      cov_names <- colnames(predictors)[-cols_to_remove]
    } else {
      X_mat <- matrix(nrow = n, ncol = 0)
      cov_names <- character(0)
    }
  } else {
    group_only <- hardhat::mold(
      stats::reformulate(group_col),
      original_data
    )
    group_predictor_names <- colnames(group_only$predictors)

    keep_cols <- !colnames(predictors) %in% group_predictor_names
    keep_cols <- keep_cols & colnames(predictors) != group_col
    X_mat <- as.matrix(predictors[, keep_cols, drop = FALSE])

    cov_names <- if (ncol(X_mat) > 0) colnames(X_mat) else character(0)
  }

  if (nrow(X_mat) != n) {
    stop(
      "Covariate matrix has ",
      nrow(X_mat),
      " rows but outcome has ",
      n,
      " rows.",
      call. = FALSE
    )
  }

  if (any(is.na(X_mat))) {
    stop(
      "Covariate matrix contains missing values. ",
      "Please handle missingness before model fitting.",
      call. = FALSE
    )
  }

  if (!all(is.finite(X_mat))) {
    stop("Covariate matrix contains non-finite values.", call. = FALSE)
  }

  n_covariates <- ncol(X_mat)
  checkmate::assert_int(n_covariates, lower = 0)

  # prepare cpp data
  cpp_data <- list(
    time = time_vec,
    event = event_vec,
    group = group_vec,
    X = X_mat,
    n = n,
    n_groups = n_groups,
    n_covariates = n_covariates
  )

  # run mcmc
  cpp_out <- cpp_exnex_gibbs(
    time = cpp_data$time,
    event = cpp_data$event,
    group = cpp_data$group,
    X = cpp_data$X,
    priors = priors,
    iter = iter,
    warmup = warmup,
    chains = chains
  )

  if (!is.data.frame(cpp_out$draws)) {
    cpp_out$draws <- as.data.frame(cpp_out$draws)
  }

  checkmate::assert_data_frame(cpp_out$draws, min.rows = 1, min.cols = 1)

  # store metadata
  clean_data <- list(
    time = time_vec,
    event = event_vec,
    group = group_vec,
    X = X_mat,
    n = n,
    n_groups = n_groups,
    n_covariates = n_covariates,
    cov_names = cov_names
  )

  new_exnex_surv(
    draws = cpp_out$draws,
    data = clean_data,
    priors = priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    blueprint = processed$blueprint
  )
}

# helpers

#' Extract time vector from outcomes
#' @keywords internal
.extract_time_vector <- function(outcomes) {
  if (inherits(outcomes[[1]], "Surv") || is.matrix(outcomes[[1]])) {
    as.numeric(outcomes[[1]][, 1])
  } else if (is.data.frame(outcomes) && ncol(outcomes) >= 2) {
    as.numeric(outcomes[[1]])
  } else {
    stop(
      "Outcome must contain time and event status. ",
      "Use `Surv(time, event)` on the left-hand side of the formula.",
      call. = FALSE
    )
  }
}

#' Extract event vector from outcomes
#' @keywords internal
.extract_event_vector <- function(outcomes) {
  if (inherits(outcomes[[1]], "Surv") || is.matrix(outcomes[[1]])) {
    as.numeric(outcomes[[1]][, 2])
  } else if (is.data.frame(outcomes) && ncol(outcomes) >= 2) {
    as.numeric(outcomes[[2]])
  } else {
    stop(
      "Outcome must contain time and event status. ",
      "Use `Surv(time, event)` on the left-hand side of the formula.",
      call. = FALSE
    )
  }
}

#' Extract group vector from one-hot encoded predictors
#' Assumes first set of consecutive columns with same prefix form the group
#' Validates proper one-hot encoding: exactly one 1 per row
#' @keywords internal
.extract_group_from_encoded <- function(predictors) {
  col_names <- colnames(predictors)
  first_name <- col_names[1]
  base_name <- gsub("[0-9]+$", "", first_name)
  group_cols <- grep(paste0("^", base_name), col_names)

  if (length(group_cols) == 0) {
    stop("Could not identify group columns in predictors.", call. = FALSE)
  }

  group_data <- as.matrix(predictors[, group_cols, drop = FALSE])
  row_sums <- rowSums(group_data)

  if (!all(row_sums == 1)) {
    bad_rows <- which(row_sums != 1)
    stop(
      "Invalid one-hot encoding: ",
      length(bad_rows),
      " row(s) do not have exactly one 1. ",
      "First bad row: ",
      bad_rows[1],
      " (sum=",
      row_sums[bad_rows[1]],
      ")",
      call. = FALSE
    )
  }

  group_vec <- apply(group_data, 1, function(row) {
    which(row == 1)
  })

  if (!is.numeric(group_vec)) {
    stop(
      "Failed to extract group vector. Result is ",
      class(group_vec)[1],
      " instead of numeric.",
      call. = FALSE
    )
  }

  return(as.numeric(group_vec))
}

#' Count number of group columns in one-hot encoded predictors
#' Uses the same validation as extract_group_from_encoded to be consistent
#' @keywords internal
.count_group_columns <- function(predictors) {
  col_names <- colnames(predictors)
  first_name <- col_names[1]
  base_name <- gsub("[0-9]+$", "", first_name)

  group_cols <- grep(paste0("^", base_name), col_names)

  if (length(group_cols) == 0) {
    stop("Could not identify group columns in predictors.", call. = FALSE)
  }

  group_data <- as.matrix(predictors[, group_cols, drop = FALSE])

  if (!all(group_data %in% c(0, 1))) {
    stop(
      "Group columns do not contain only 0/1 values (not one-hot encoded).",
      call. = FALSE
    )
  }

  return(length(group_cols))
}
