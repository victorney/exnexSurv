#' Bridge connecting hardhat processed data to the C++ Gibbs Sampler
#'
#' Orchestrates data preparation, validation, and passes it to the C++ MCMC engine.
#'
#' The group variable (specified by group_col) is identified from the original
#' data. Additional variables are treated as covariates in the linear predictor.
#'
#' @param processed A list produced by `hardhat::mold()`.
#' @param pooling Pooling mode; see \code{\link{pooling_surv}}().
#' @param priors Optional named list of prior hyperparameters; see
#'   \code{pooling_surv()}.
#' @param iter Total number of MCMC iterations.
#' @param warmup Number of warmup iterations.
#' @param thin Thinning interval for the post-warmup draws.
#' @param chains Number of chains to run.
#' @param parallel_chains Number of chains to run in parallel at the R level.
#' @param group_col Name of the original group column.
#' @param original_data The original data frame (before hardhat processing).
#' @param seed Optional seed used to deterministically derive one seed per chain.
#'
#' @return An `pooling_surv` object.
#' @keywords internal
pooling_surv_bridge <- function(
  processed,
  pooling,
  priors,
  iter,
  warmup,
  thin,
  chains,
  parallel_chains,
  group_col,
  original_data,
  seed = NULL,
  verbose = TRUE
) {
  predictors <- processed$predictors
  outcomes <- processed$outcomes

  # extract outcomes
  time_vec <- .extract_time_vector(outcomes)
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

  if (any(is.na(group_vec))) {
    stop("Group variable contains missing values.", call. = FALSE)
  }

  n_groups <- max(group_vec)
  checkmate::assert_int(n_groups, lower = 1)

  # extract covariates
  group_only <- hardhat::mold(
    stats::reformulate(group_col),
    original_data
  )
  group_predictor_names <- colnames(group_only$predictors)

  keep_cols <- !colnames(predictors) %in% group_predictor_names
  keep_cols <- keep_cols & colnames(predictors) != group_col
  X_mat <- as.matrix(predictors[, keep_cols, drop = FALSE])

  cov_names <- if (ncol(X_mat) > 0) colnames(X_mat) else character(0)

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

  # Warn once about priors that are irrelevant for the selected pooling mode.
  .warn_ignored_priors(priors, pooling)

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

  chain_seeds <- .make_chain_seeds(chains = chains, seed = seed)

  chain_runs <- if (chains == 1 || !parallel_chains) {
    run_chains <- function(progressor = NULL) {
      lapply(
        seq_len(chains),
        function(chain_id) {
          hook <- if (!is.null(progressor)) {
            function(msg) progressor(msg)
          } else {
            NULL
          }
          .run_single_chain_pooling(
            cpp_data = cpp_data,
            priors = priors,
            pooling = pooling,
            iter = iter,
            warmup = warmup,
            thin = thin,
            seed = chain_seeds[chain_id],
            chain_label = as.character(chain_id),
            verbose = FALSE,
            progress_hook = hook
          )
        }
      )
    }
    if (verbose) {
      progressr::with_progress(
        run_chains(progressr::progressor(steps = 10L * chains)),
        enable = TRUE
      )
    } else {
      run_chains()
    }
  } else {
    .run_chains_parallel_pooling(
      cpp_data = cpp_data,
      priors = priors,
      pooling = pooling,
      verbose = verbose,
      iter = iter,
      warmup = warmup,
      thin = thin,
      chains = chains,
      chain_seeds = chain_seeds
    )
  }

  chain_draws <- lapply(chain_runs, function(x) x$draws)
  chain_diagnostics <- lapply(chain_runs, function(x) x$diagnostics)

  draws <- do.call(rbind, chain_draws)
  rownames(draws) <- NULL

  if (!is.data.frame(draws)) {
    draws <- as.data.frame(draws)
  }

  checkmate::assert_data_frame(draws, min.rows = 1, min.cols = 1)

  # store metadata
  clean_data <- list(
    time = time_vec,
    event = event_vec,
    group = group_vec,
    X = X_mat,
    n = n,
    n_groups = n_groups,
    n_covariates = n_covariates,
    cov_names = cov_names,
    chain_seeds = chain_seeds
  )

  new_pooling_surv(
    draws = draws,
    data = clean_data,
    pooling = pooling,
    priors = priors,
    resolved_priors = chain_diagnostics[[1L]]$resolved_priors,
    iter = iter,
    warmup = warmup,
    thin = thin,
    chains = chains,
    blueprint = processed$blueprint
  )
}

# helpers

#' Warn about priors that are ignored in the selected pooling mode
#' @keywords internal
.warn_ignored_priors <- function(priors, pooling) {
  used <- switch(pooling,
    exnex = c(
      "a_sigma", "b_sigma", "a_tau", "b_tau", "p_mix",
      "m_mu", "v_mu", "m_nex", "v_nex", "v_beta"
    ),
    complete = c("a_sigma", "b_sigma", "m_mu", "v_mu", "v_beta"),
    none = c("a_sigma", "b_sigma", "m_nex", "v_nex", "v_beta")
  )

  supplied <- names(priors)
  ignored <- supplied[!supplied %in% used]

  if (length(ignored) > 0) {
    warning(
      "Priors ", paste0("'", ignored, "'", collapse = ", "),
      " were specified but are ignored with pooling = \"", pooling, "\".",
      call. = FALSE
    )
  }
}

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

#' Derive one deterministic seed per chain (restores the global RNG state)
#' @keywords internal
.make_chain_seeds <- function(chains, seed = NULL) {
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit(
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      },
      add = TRUE
    )
    set.seed(seed)
  }
  sample.int(.Machine$integer.max, size = chains, replace = FALSE)
}

#' Signal progress without perturbing the worker RNG stream
#'
#' `progressor()` (and condition signalling in general) draws random
#' condition identifiers inside the worker; restoring `.Random.seed` around
#' the call keeps the MCMC draws independent of the progress reporting.
#' @keywords internal
.wrap_progress_hook <- function(hook) {
  function(msg) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    hook(msg)
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv, inherits = FALSE)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    invisible(NULL)
  }
}

#' Run a single chain by calling the C++ kernel once
#'
#' `rng_kind` pins the RNG configuration (future workers may default to a
#' different normal-kind sampler, which would make the same seed produce
#' different draws than sequential execution).
#' @keywords internal
.run_single_chain_pooling <- function(cpp_data, priors, pooling, iter, warmup, seed, thin = 1L, verbose = TRUE, chain_label = "", progress_hook = NULL, rng_kind = NULL) {
  # future workers may default to a different normal-kind sampler than the
  # master session; pin the master's RNG kind so the same seed produces
  # identical draws in sequential and parallel execution
  if (!is.null(rng_kind)) {
    do.call(RNGkind, as.list(rng_kind))
  }
  set.seed(seed)

  out <- cpp_exnex_gibbs(
    time = cpp_data$time,
    event = cpp_data$event,
    group = cpp_data$group,
    X = cpp_data$X,
    priors = priors,
    pooling = pooling,
    verbose = verbose,
    iter = iter,
    warmup = warmup,
    thin = thin,
    chains = 1L,
    chain_label = chain_label,
    progress_hook = if (is.null(progress_hook)) NULL else .wrap_progress_hook(progress_hook)
  )

  list(
    draws = as.data.frame(out$draws),
    diagnostics = out$diagnostics
  )
}

#' Run all chains concurrently on future workers
#' @keywords internal
.run_chains_parallel_pooling <- function(
  cpp_data,
  priors,
  pooling,
  verbose,
  iter,
  warmup,
  thin,
  chains,
  chain_seeds
) {
  old_plan <- future::plan()
  old_rng_misuse <- getOption("future.rng.onMisuse")
  master_rng_kind <- RNGkind()
  future::plan(future::multisession, workers = chains)
  options(future.rng.onMisuse = "ignore")
  on.exit({
    future::plan(old_plan)
    options(future.rng.onMisuse = old_rng_misuse)
  }, add = TRUE)

  relay_chains <- function(progressor = NULL) {
    future.apply::future_lapply(
      X = seq_len(chains),
      FUN = function(chain_id, cpp_data, priors, pooling, iter, warmup, thin, chain_seeds, master_rng_kind, progressor) {
        hook <- if (!is.null(progressor)) {
          function(msg) progressor(msg)
        } else {
          NULL
        }
        .run_single_chain_pooling(
          cpp_data = cpp_data,
          priors = priors,
          pooling = pooling,
          iter = iter,
          warmup = warmup,
          thin = thin,
          seed = chain_seeds[chain_id],
          verbose = FALSE,
          chain_label = as.character(chain_id),
          progress_hook = hook,
          rng_kind = master_rng_kind
        )
      },
      cpp_data = cpp_data,
      priors = priors,
      pooling = pooling,
      iter = iter,
      warmup = warmup,
      thin = thin,
      chain_seeds = chain_seeds,
      master_rng_kind = master_rng_kind,
      progressor = progressor
    )
  }

  # Workers cannot write to the master console; progress is a live
  # progress bar via `progressr` (future relays worker progress conditions
  # continuously). No per-chain text lines are printed in parallel mode.
  if (verbose) {
    progressr::with_progress(
      relay_chains(progressr::progressor(steps = 10L * chains)),
      enable = TRUE
    )
  } else {
    relay_chains()
  }
}
