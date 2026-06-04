#' Main Gibbs Sampler for EXNEX Survival Models
#'
#' Internal R wrapper around the compiled C++ sampler.
#'
#' This function is intentionally kept in a regular R file so the public
#' documentation and the R-side calling contract remain separate from the
#' generated Rcpp bootstrap code.
#'
#' @param time Vector of observed follow-up times (n-vector, all > 0)
#' @param event Vector of event indicators (n-vector, 0 or 1)
#' @param group Vector of group assignments (n-vector, integers 1 to K)
#' @param X Matrix of covariates (n x P). Can be empty (n x 0) if no covariates.
#' @param priors List with prior specifications (placeholder for now)
#' @param iter Total number of MCMC iterations
#' @param warmup Number of iterations to discard (note: parameter is (iter - warmup))
#' @param chains Number of independent chains to run
#'
#' @return List containing posterior draws and diagnostics from the C++ layer.
#' @keywords internal
cpp_exnex_gibbs <- function(time, event, group, X, priors, iter, warmup, chains) {
  .Call(`_exnexSurv_cpp_exnex_gibbs`, time, event, group, X, priors, iter, warmup, chains)
}
