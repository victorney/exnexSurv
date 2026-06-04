#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List cpp_exnex_gibbs(const arma::vec &time, const arma::vec &event,
                           const arma::vec &group, const arma::mat &X,
                           Rcpp::List priors, const int &iter,
                           const int &warmup, const int &chains) {

  int n = time.n_elem;

  // ========== VALIDATE MCMC PARAMETERS ==========
  // Check for integer underflow: warmup >= iter causes n_samples to underflow
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

  // ========== VALIDATE VECTOR LENGTHS ==========
  // Ensure all vectors have matching length
  if (event.n_elem != n) {
    Rcpp::stop("event must have the same length as time. Got time.n_elem=" +
               std::to_string(n) + ", event.n_elem=" + std::to_string(event.n_elem));
  }
  if (group.n_elem != n) {
    Rcpp::stop("group must have the same length as time. Got time.n_elem=" +
               std::to_string(n) + ", group.n_elem=" + std::to_string(group.n_elem));
  }
  if (X.n_rows > 0 && X.n_rows != n) {
    Rcpp::stop("X must have the same number of rows as time. Got time.n_elem=" +
               std::to_string(n) + ", X.n_rows=" + std::to_string(X.n_rows));
  }

  // ========== VALIDATE GROUP VECTOR ==========
  // Group must contain integers >= 1
  double min_group = arma::min(group);
  double max_group = arma::max(group);
  
  if (min_group < 1.0) {
    Rcpp::stop("All group values must be >= 1. Found minimum: " +
               std::to_string(min_group));
  }
  
  // Check that group values are actually integers
  for (int i = 0; i < group.n_elem; ++i) {
    if (group(i) != std::floor(group(i))) {
      Rcpp::stop("Group must contain integer values only. Found non-integer: " +
                 std::to_string(group(i)) + " at position " + std::to_string(i + 1));
    }
  }

  // ========== EXTRACT DIMENSIONS (NOW SAFE) ==========
  int K = static_cast<int>(max_group); // Number of groups
  int P = X.n_cols; // Number of covariates (can be 0)
  int n_samples = iter - warmup; // Number of samples to retain
  int n_cols_out = K + P + 1; // Total columns: theta + beta + sigma

  if (K < 1) {
    Rcpp::stop("At least one group must be present.");
  }
  
  if (n_samples <= 0) {
    Rcpp::stop("No samples to retain after warmup. n_samples=" + 
               std::to_string(n_samples));
  }

  // each row is an MCMC iteration (after warmup)
  // columns: theta_1, ..., theta_K, beta_1, ..., beta_P, sigma
  arma::mat theta_samples(n_samples, K, arma::fill::zeros);
  arma::mat beta_samples(n_samples, P, arma::fill::zeros);
  arma::vec sigma_samples(n_samples, arma::fill::ones);

  arma::mat raw_draws(n_samples, n_cols_out, arma::fill::zeros);

  // MCMC placeholder
  //...

  // assign theta columns (0 to K-1)
  raw_draws.cols(0, K - 1) = theta_samples;

  // assign beta columns (K to K+P-1), if P > 0
  if (P > 0) {
    raw_draws.cols(K, K + P - 1) = beta_samples;
  }

  // assign sigma column (K+P)
  raw_draws.col(K + P) = sigma_samples;

  Rcpp::NumericMatrix draws_out = Rcpp::wrap(raw_draws);
  Rcpp::CharacterVector col_names(n_cols_out);

  int col_idx = 0;

  // theta_1, theta_2, ..., theta_K
  for (int k = 0; k < K; ++k) {
    col_names[col_idx++] = "theta_" + std::to_string(k + 1);
  }

  // beta_1, beta_2, ..., beta_P (if P > 0)
  for (int p = 0; p < P; ++p) {
    col_names[col_idx++] = "beta_" + std::to_string(p + 1);
  }

  // sigma
  col_names[col_idx] = "sigma";
  Rcpp::colnames(draws_out) = col_names;

  // prepare diagnostics list
  Rcpp::List diagnostics =
      Rcpp::List::create(Rcpp::Named("n_obs") = n, Rcpp::Named("n_groups") = K,
                         Rcpp::Named("n_covariates") = P,
                         Rcpp::Named("n_censored") = arma::sum(1 - event),
                         Rcpp::Named("n_observed") = arma::sum(event));

  // return results
  return Rcpp::List::create(
      Rcpp::Named("draws") = draws_out, Rcpp::Named("priors") = priors,
      Rcpp::Named("iter") = iter, Rcpp::Named("warmup") = warmup,
      Rcpp::Named("chains") = chains, Rcpp::Named("diagnostics") = diagnostics);
}
