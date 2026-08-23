#' @keywords internal
#' @description
#' Fits Bayesian Exchangeable--Non-Exchangeable (EXNEX) survival models for
#' right-censored log-normal data in basket trials using a fast,
#' data-augmented Gibbs sampler written in C++ (Rcpp/RcppArmadillo).
#'
#' For patient \eqn{i}, \eqn{\log T_i=\theta_{g[i]}+X_i^{\mathsf T}\beta+
#' \varepsilon_i} with \eqn{\varepsilon_i\sim\mathcal N(0,\sigma^2)}. A
#' latent indicator \eqn{Z_j} selects, for each basket \eqn{j}, whether its
#' effect is drawn from the exchangeable component
#' \eqn{\mathcal N(\mu,\tau^2)} or from a basket-specific non-exchangeable
#' prior \eqn{\mathcal N(m_{0j},v_{0j})}. Censored event times are handled by
#' data augmentation (imputing truncated-Normal latent log-times), which makes
#' every full conditional conjugate and the systematic Gibbs scan exact.
#'
#' The main entry point is \code{\link{exnex_surv}()}. See the package
#' vignette \emph{The EXNEX Model, Priors, and Data Augmentation} for a full
#' presentation of the model, the priors, and the sampler.
#'
#' @references
#' Neuenschwander, B., Wandel, S., Roychoudhury, S., & Bailey, S. (2016).
#' Robust exchangeability designs for early phase clinical trials with multiple
#' strata. \emph{Pharmaceutical Statistics}, 15(2), 123--134.
#'
#' Tanner, M. A., & Wong, W. H. (1987). The calculation of posterior
#' distributions by data augmentation. \emph{Journal of the American
#' Statistical Association}, 82(398), 528--540.
#'
#' @seealso \code{\link{exnex_surv}}, \code{\link{summary.exnex_surv}},
#'   \code{\link{print.exnex_surv}}, \code{\link{plot.exnex_surv}}
"_PACKAGE"
#' @useDynLib exnexSurv, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats pnorm rnorm runif
NULL
