#include <RcppArmadillo.h>
#include <cmath>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]

//' Main Gibbs Sampler for EXNEX Survival Models
//'
//' Implements MCMC sampling for Bayesian EXNEX model with right-censored
//' log-normal survival data. Supports optional covariates.
//'
//' @param time Vector of observed follow-up times (n-vector, all > 0)
//' @param event Vector of event indicators (n-vector, 0 or 1)
//' @param group Vector of group assignments (n-vector, integers 1 to K)
//' @param X Matrix of covariates (n x P). Can be empty (n x 0) if no covariates.
//' @param priors List with prior specifications
//' @param iter Total number of MCMC iterations
//' @param warmup Number of iterations to discard
//' @param chains Number of independent chains to run
//'
//' @return List containing posterior draws, priors, metadata, and diagnostics
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_exnex_gibbs(const arma::vec &time, const arma::vec &event,
                           const arma::vec &group, const arma::mat &X,
                           Rcpp::List priors, const int &iter,
                           const int &warmup, const int &chains) {

  int n = time.n_elem;

  if (n == 0) {
    Rcpp::stop("time must have positive length.");
  }

  if (!time.is_finite() || arma::any(time <= 0)) {
    Rcpp::stop("All survival times must be positive and finite.");
  }
  // validate mcmc parameters
  if (iter <= 0) {
    Rcpp::stop("iter must be positive.");
  }
  if (warmup < 0) {
    Rcpp::stop("warmup must be non-negative.");
  }
  if (warmup >= iter) {
    Rcpp::stop("warmup must be strictly less than iter. Got iter=" +
               std::to_string(iter) + ", warmup=" + std::to_string(warmup));
  }

  // validate vector lengths
  if (event.n_elem != n) {
    Rcpp::stop("event must have the same length as time. Got time.n_elem=" +
               std::to_string(n) + ", event.n_elem=" + std::to_string(event.n_elem));
  }
  if (arma::any((event != 0.0) % (event != 1.0))) {
    Rcpp::stop("event must contain only 0/1 values.");
  }
  if (group.n_elem != n) {
    Rcpp::stop("group must have the same length as time. Got time.n_elem=" +
               std::to_string(n) + ", group.n_elem=" + std::to_string(group.n_elem));
  }
  if (X.n_rows != n) {
    Rcpp::stop("X must have the same number of rows as time. Got time.n_elem=" +
               std::to_string(n) + ", X.n_rows=" + std::to_string(X.n_rows));
  }

  // validate group vector
  double min_group = arma::min(group);
  double max_group = arma::max(group);
  
  if (min_group < 1.0) {
    Rcpp::stop("All group values must be >= 1. Found minimum: " +
               std::to_string(min_group));
  }
  
  for (int i = 0; i < group.n_elem; ++i) {
    if (group(i) != std::floor(group(i))) {
      Rcpp::stop("Group must contain integer values only. Found non-integer: " +
                 std::to_string(group(i)) + " at position " + std::to_string(i + 1));
    }
  }

  // extract dimensions
  int K = static_cast<int>(max_group); 
  int P = X.n_cols; 
  int n_samples = iter - warmup; 
  int n_cols_out = K + P + 1; // K + P + sigma

  if (K < 1) {
    Rcpp::stop("At least one group must be present.");
  }
  
  if (n_samples <= 0) {
    Rcpp::stop("No samples to retain after warmup. n_samples=" + 
               std::to_string(n_samples));
  }

  arma::mat theta_samples(n_samples, K, arma::fill::zeros);
  arma::mat beta_samples(n_samples, P, arma::fill::zeros);
  arma::vec sigma_samples(n_samples, arma::fill::ones);

  arma::mat raw_draws(n_samples, n_cols_out, arma::fill::zeros);

  // mcmc placeholder
  //...

  raw_draws.cols(0, K - 1) = theta_samples;

  if (P > 0) {
    raw_draws.cols(K, K + P - 1) = beta_samples;
  }

  raw_draws.col(K + P) = sigma_samples;

  Rcpp::NumericMatrix draws_out = Rcpp::wrap(raw_draws);
  Rcpp::CharacterVector col_names(n_cols_out);

  int col_idx = 0;

  for (int k = 0; k < K; ++k) {
    col_names[col_idx++] = "theta_" + std::to_string(k + 1);
  }

  for (int p = 0; p < P; ++p) {
    col_names[col_idx++] = "beta_" + std::to_string(p + 1);
  }

  col_names[col_idx] = "sigma";
  Rcpp::colnames(draws_out) = col_names;

  Rcpp::List diagnostics =
      Rcpp::List::create(Rcpp::Named("n_obs") = n, Rcpp::Named("n_groups") = K,
                         Rcpp::Named("n_covariates") = P,
                         Rcpp::Named("n_censored") = arma::sum(1 - event),
                         Rcpp::Named("n_observed") = arma::sum(event));

  return Rcpp::List::create(
      Rcpp::Named("draws") = draws_out, Rcpp::Named("priors") = priors,
      Rcpp::Named("iter") = iter, Rcpp::Named("warmup") = warmup,
      Rcpp::Named("chains") = chains, Rcpp::Named("diagnostics") = diagnostics);
}
