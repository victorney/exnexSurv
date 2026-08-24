#include <RcppArmadillo.h>
#include <cmath>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]

double draw_inverse_gamma(double shape, double rate) {
  if (shape <= 0.0 || rate <= 0.0) {
    Rcpp::stop("Inverse-gamma parameters must be positive.");
  }

  return 1.0 / R::rgamma(shape, 1.0 / rate);
}

double draw_left_trunc_normal(double lower, double mean, double sd) {
  // Inverse-CDF sampler for a normal truncated below at `lower`.
  double eps = 1e-12;

  if (!std::isfinite(lower) || !std::isfinite(mean) || !std::isfinite(sd) || sd <= 0.0) {
    Rcpp::stop("Truncated normal parameters must be finite with sd > 0.");
  }

  double z_lower = (lower - mean) / sd;
  double log_tail = R::pnorm5(z_lower, 0.0, 1.0, 0, 1);

  if (!std::isfinite(log_tail)) {
    return lower + eps;
  }

  double u = R::runif(0.0, 1.0);
  if (u < eps) {
    u = eps;
  }

  double z = R::qnorm5(std::log(u) + log_tail, 0.0, 1.0, 0, 1);
  double out = mean + sd * z;

  if (!std::isfinite(out) || out < lower) {
    out = lower + eps;
  }

  return out;
}

arma::vec draw_normal_from_precision(const arma::mat& precision, const arma::vec& rhs) {
  // Draw from N(Q^{-1} rhs, Q^{-1}) using the precision matrix Q.
  if (rhs.n_elem == 0) {
    return arma::vec();
  }

  arma::mat precision_use = 0.5 * (precision + precision.t());
  arma::mat chol_upper;
  double jitter = 1e-8;
  bool ok = false;

  for (int attempt = 0; attempt < 8; ++attempt) {
    ok = arma::chol(chol_upper, precision_use);
    if (ok) {
      break;
    }

    precision_use.diag() += jitter;
    jitter *= 10.0;
  }

  if (!ok) {
    Rcpp::stop("Failed to compute a stable Cholesky decomposition.");
  }

  arma::vec tmp = arma::solve(arma::trimatl(chol_upper.t()), rhs);
  arma::vec mean = arma::solve(arma::trimatu(chol_upper), tmp);
  arma::vec z = arma::zeros<arma::vec>(rhs.n_elem);

  for (arma::uword i = 0; i < z.n_elem; ++i) {
    z(i) = R::rnorm(0.0, 1.0);
  }

  arma::vec noise = arma::solve(arma::trimatu(chol_upper), z);
  return mean + noise;
}

//' Main Gibbs Sampler for EXNEX Survival Models
//'
//' Implements Gibbs sampling for Bayesian EXNEX right-censored
//' log-normal survival models with optional covariates.
//'
//' @param time Vector of observed follow-up times (n-vector, all > 0)
//' @param event Vector of event indicators (n-vector, 0 or 1)
//' @param group Vector of group assignments (n-vector, integers 1 to K)
//' @param X Matrix of covariates (n x P). Can be empty (n x 0) if no covariates.
//' @param priors Optional named list of prior hyperparameters. Supported
//'   fields: \code{a_sigma}, \code{b_sigma}, \code{a_tau}, \code{b_tau}
//'   (inverse-Gamma shape and scale), \code{p_mix} (EXNEX mixture weight),
//'   \code{m_mu}, \code{v_mu} (exchangeable mean prior), \code{m_nex},
//'   \code{v_nex} (nonexchangeable component), \code{v_beta} (variance of the
//'   regression-coefficient prior). \code{p_mix}, \code{m_nex}, and
//'   \code{v_nex} each accept either a scalar, replicated across baskets, or a
//'   numeric vector of length K with one value per basket, matching the
//'   basket-specific notation \eqn{p_{\mathrm{exch},j}}, \eqn{m_{0j}},
//'   \eqn{v_{0j}} of the model. Absent fields keep the defaults; unknown
//'   fields are ignored.
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
  Rcpp::RNGScope scope;

  double a_sigma = 2.0;
  double b_sigma = 2.0;
  double a_tau = 2.0;
  double b_tau = 2.0;
  double p_mix = 0.5;
  double m_mu = 0.0;
  double v_mu = 1e4;
  double m_nex = 0.0;
  double v_nex = 1e4;
  double v_beta = 1e4;
  double init_censor_offset = 0.1;

  // Track whether the basket-specific priors were supplied as vectors so that
  // the resolved-priors report preserves the user's input shape.
  bool p_mix_vector = false;
  bool m_nex_vector = false;
  bool v_nex_vector = false;

  // Read optional prior overrides from the R list. Absent fields keep the
  // defaults above; unknown fields are ignored for backward compatibility.
  auto read_prior = [](const Rcpp::List& priors, const char* name,
                       double default_value) {
    if (!priors.containsElementNamed(name)) {
      return default_value;
    }
    SEXP element = priors[name];
    if (Rf_isNull(element) || !Rf_isNumeric(element) || Rf_length(element) != 1) {
      Rcpp::stop(std::string("Prior '") + name + "' must be a single numeric value.");
    }
    double value = Rcpp::as<double>(element);
    if (!std::isfinite(value)) {
      Rcpp::stop(std::string("Prior '") + name + "' must be finite.");
    }
    return value;
  };

  a_sigma = read_prior(priors, "a_sigma", a_sigma);
  b_sigma = read_prior(priors, "b_sigma", b_sigma);
  a_tau = read_prior(priors, "a_tau", a_tau);
  b_tau = read_prior(priors, "b_tau", b_tau);
  // p_mix, m_nex and v_nex may be scalar or per-basket vectors; the scalar
  // defaults are read below by read_prior_vec, which replicates them.
  m_mu = read_prior(priors, "m_mu", m_mu);
  v_mu = read_prior(priors, "v_mu", v_mu);
  v_beta = read_prior(priors, "v_beta", v_beta);

  if (a_sigma <= 0.0 || b_sigma <= 0.0 || a_tau <= 0.0 || b_tau <= 0.0) {
    Rcpp::stop("Inverse-gamma prior parameters must be positive.");
  }
  if (v_mu <= 0.0 || v_beta <= 0.0) {
    Rcpp::stop("Prior variances must be positive.");
  }

  int n = time.n_elem;

  if (n == 0) {
    Rcpp::stop("time must have positive length.");
  }

  if (!time.is_finite() || arma::any(time <= 0)) {
    Rcpp::stop("All survival times must be positive and finite.");
  }

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

  if (event.n_elem != static_cast<arma::uword>(n)) {
    Rcpp::stop("event must have the same length as time. Got time.n_elem=" +
               std::to_string(n) + ", event.n_elem=" + std::to_string(event.n_elem));
  }

  if (arma::any((event != 0.0) % (event != 1.0))) {
    Rcpp::stop("event must contain only 0/1 values.");
  }

  if (group.n_elem != static_cast<arma::uword>(n)) {
    Rcpp::stop("group must have the same length as time. Got time.n_elem=" +
               std::to_string(n) + ", group.n_elem=" + std::to_string(group.n_elem));
  }

  if (X.n_rows != static_cast<arma::uword>(n)) {
    Rcpp::stop("X must have the same number of rows as time. Got time.n_elem=" +
               std::to_string(n) + ", X.n_rows=" + std::to_string(X.n_rows));
  }

  double min_group = arma::min(group);
  double max_group = arma::max(group);

  if (min_group < 1.0) {
    Rcpp::stop("All group values must be >= 1. Found minimum: " +
               std::to_string(min_group));
  }

  for (arma::uword i = 0; i < group.n_elem; ++i) {
    if (group(i) != std::floor(group(i))) {
      Rcpp::stop("Group must contain integer values only. Found non-integer: " +
                 std::to_string(group(i)) + " at position " + std::to_string(i + 1));
    }
  }

  int K = static_cast<int>(max_group);
  int P = X.n_cols;
  int n_samples = iter - warmup;
  int n_cols_out = K + P + 1;

  // Read basket-specific priors: a scalar is replicated across the K baskets;
  // a numeric vector of length K assigns one value per basket. The scalars
  // read above serve as defaults and as the replicated value.
  auto read_prior_vec = [&](const char* name, double default_value,
                            bool& was_vector) {
    arma::vec out(K);
    if (!priors.containsElementNamed(name)) {
      out.fill(default_value);
      return out;
    }
    SEXP element = priors[name];
    if (Rf_isNull(element) || !Rf_isNumeric(element)) {
      Rcpp::stop(std::string("Prior '") + name +
                 "' must be a single numeric value or a numeric vector of "
                 "length K (one value per basket).");
    }
    if (Rf_length(element) == 1) {
      double value = Rcpp::as<double>(element);
      if (!std::isfinite(value)) {
        Rcpp::stop(std::string("Prior '") + name + "' must be finite.");
      }
      out.fill(value);
      return out;
    }
    if (Rf_length(element) == K) {
      out = Rcpp::as<arma::vec>(element);
      if (!out.is_finite()) {
        Rcpp::stop(std::string("Prior '") + name + "' must be finite.");
      }
      was_vector = true;
      return out;
    }
    Rcpp::stop(std::string("Prior '") + name + "' must be a single numeric "
               "value or a numeric vector of length K = " +
               std::to_string(K) + "; got length " +
               std::to_string(Rf_length(element)) + ".");
    return out;
  };

  arma::vec p_mix_vec = read_prior_vec("p_mix", p_mix, p_mix_vector);
  arma::vec m_nex_vec = read_prior_vec("m_nex", m_nex, m_nex_vector);
  arma::vec v_nex_vec = read_prior_vec("v_nex", v_nex, v_nex_vector);

  // Keep the scalar mirrors in sync so that the resolved-priors report
  // reflects scalar input; vector input is reported as a vector.
  if (!p_mix_vector) p_mix = p_mix_vec(0);
  if (!m_nex_vector) m_nex = m_nex_vec(0);
  if (!v_nex_vector) v_nex = v_nex_vec(0);

  if (arma::any(p_mix_vec <= 0.0) || arma::any(p_mix_vec >= 1.0)) {
    Rcpp::stop("Prior mixture weight 'p_mix' must lie strictly between 0 and 1 "
               "for every basket.");
  }
  if (arma::any(v_nex_vec <= 0.0)) {
    Rcpp::stop("Prior variance 'v_nex' must be positive for every basket.");
  }

  arma::vec log_time = arma::log(time);
  arma::uvec group_index = arma::conv_to<arma::uvec>::from(group - 1.0);

  // Sparse membership map used for subgroup aggregation.
  arma::sp_mat group_map(n, K);
  for (int i = 0; i < n; ++i) {
    group_map(i, group_index(i)) = 1.0;
  }

  arma::vec group_counts = group_map.t() * arma::ones<arma::vec>(n);

  arma::vec y_star = log_time;
  arma::uvec censored_index = arma::find(event == 0.0);

  // Start censored values slightly above the censoring threshold.
  if (!censored_index.is_empty()) {
    y_star.elem(censored_index) = log_time.elem(censored_index) + init_censor_offset;
  }

  arma::vec beta = arma::zeros<arma::vec>(P);
  arma::vec theta = arma::zeros<arma::vec>(K);
  arma::ivec z = arma::ones<arma::ivec>(K);

  double mu = arma::mean(y_star);
  for (int j = 0; j < K; ++j) {
    arma::uvec idx = arma::find(group_index == static_cast<arma::uword>(j));
    if (!idx.is_empty()) {
      theta(j) = arma::mean(y_star.elem(idx));
    } else {
      theta(j) = mu;
    }
  }

  double sigma2 = arma::var(y_star);
  if (!std::isfinite(sigma2) || sigma2 <= 0.0) {
    sigma2 = 1.0;
  }

  double tau2 = arma::var(theta);
  if (!std::isfinite(tau2) || tau2 <= 0.0) {
    tau2 = 1.0;
  }

  arma::mat draws = arma::zeros<arma::mat>(n_samples, n_cols_out);
  arma::mat XtX;
  arma::mat beta_prior_precision;

  if (P > 0) {
    XtX = X.t() * X;
    beta_prior_precision = arma::eye(P, P) / v_beta;
  }

  for (int iter_idx = 0; iter_idx < iter; ++iter_idx) {
    // Current linear predictor from fixed effects.
    arma::vec xb = arma::zeros<arma::vec>(n);
    if (P > 0) {
      xb = X * beta;
    }

    // Impute latent log-times for censored observations.
    for (arma::uword h = 0; h < censored_index.n_elem; ++h) {
      arma::uword i = censored_index(h);
      double mean_i = theta(group_index(i)) + xb(i);
      double sd_i = std::sqrt(sigma2);
      y_star(i) = draw_left_trunc_normal(log_time(i), mean_i, sd_i);
    }

    arma::vec theta_by_obs = theta.elem(group_index);
    arma::vec resid = y_star - theta_by_obs - xb;

    // Update residual variance.
    double shape_sigma = a_sigma + 0.5 * n;
    double rate_sigma = b_sigma + 0.5 * arma::dot(resid, resid);
    sigma2 = draw_inverse_gamma(shape_sigma, rate_sigma);

    if (P > 0) {
      // Skip this block entirely when there are no covariates.
      arma::vec y_tilde = y_star - theta_by_obs;
      arma::mat beta_precision = XtX / sigma2 + beta_prior_precision;
      arma::vec beta_rhs = X.t() * y_tilde / sigma2;
      beta = draw_normal_from_precision(beta_precision, beta_rhs);
      xb = X * beta;
    } else {
      xb.zeros();
    }

    arma::vec subgroup_resid = y_star - xb;

    // Sum residual contributions within each subgroup.
    arma::vec sum_R = group_map.t() * subgroup_resid;

    for (int j = 0; j < K; ++j) {
      double prior_mean = m_nex_vec(j);
      double prior_var = v_nex_vec(j);

      if (z(j) == 1) {
        prior_mean = mu;
        prior_var = tau2;
      }

      double precision = group_counts(j) / sigma2 + 1.0 / prior_var;
      double variance = 1.0 / precision;
      double mean = variance * (sum_R(j) / sigma2 + prior_mean / prior_var);
      theta(j) = R::rnorm(mean, std::sqrt(variance));
    }

    for (int j = 0; j < K; ++j) {
      // Compute EX/NEX weights on the log scale for stability.
      double log_w1 = std::log(p_mix_vec(j)) + R::dnorm4(theta(j), mu, std::sqrt(tau2), 1);
      double log_w0 = std::log1p(-p_mix_vec(j)) + R::dnorm4(theta(j), m_nex_vec(j), std::sqrt(v_nex_vec(j)), 1);

      double max_log = std::max(log_w1, log_w0);
      double w1 = std::exp(log_w1 - max_log);
      double w0 = std::exp(log_w0 - max_log);
      double prob_ex = w1 / (w1 + w0);

      z(j) = R::runif(0.0, 1.0) < prob_ex ? 1 : 0;
    }

    arma::vec z_vec = arma::conv_to<arma::vec>::from(z);

    // Hyperparameter updates only use the exchangeable groups.
    double k_exchangeable = arma::accu(z_vec);
    double sum_theta_ex = arma::dot(z_vec, theta);

    double mu_precision = k_exchangeable / tau2 + 1.0 / v_mu;
    double mu_variance = 1.0 / mu_precision;
    double mu_mean = mu_variance * (sum_theta_ex / tau2 + m_mu / v_mu);
    mu = R::rnorm(mu_mean, std::sqrt(mu_variance));

    double shape_tau = a_tau + 0.5 * k_exchangeable;
    double rate_tau = b_tau + 0.5 * arma::dot(z_vec, arma::square(theta - mu));
    tau2 = draw_inverse_gamma(shape_tau, rate_tau);

    if (iter_idx >= warmup) {
      // Save only post-warmup draws.
      int draw_idx = iter_idx - warmup;
      draws(draw_idx, arma::span(0, K - 1)) = theta.t();

      if (P > 0) {
        draws(draw_idx, arma::span(K, K + P - 1)) = beta.t();
      }

      draws(draw_idx, K + P) = sigma2;
    }
  }

  Rcpp::NumericMatrix draws_out = Rcpp::wrap(draws);
  Rcpp::CharacterVector col_names(n_cols_out);
  int col_idx = 0;

  for (int k = 0; k < K; ++k) {
    col_names[col_idx] = "theta_" + std::to_string(k + 1);
    col_idx += 1;
  }

  for (int p = 0; p < P; ++p) {
    col_names[col_idx] = "beta_" + std::to_string(p + 1);
    col_idx += 1;
  }

  col_names[col_idx] = "sigma2";
  Rcpp::colnames(draws_out) = col_names;

  Rcpp::List resolved_priors = Rcpp::List::create(
    Rcpp::Named("a_sigma") = a_sigma,
    Rcpp::Named("b_sigma") = b_sigma,
    Rcpp::Named("a_tau") = a_tau,
    Rcpp::Named("b_tau") = b_tau,
    Rcpp::Named("p_mix") = p_mix,
    Rcpp::Named("m_mu") = m_mu,
    Rcpp::Named("v_mu") = v_mu,
    Rcpp::Named("m_nex") = m_nex,
    Rcpp::Named("v_nex") = v_nex,
    Rcpp::Named("v_beta") = v_beta
  );

  // Report basket-specific vector priors in their full form when supplied.
  // (arma::vec wraps with a dim attribute; rebuild as a plain R vector.)
  if (p_mix_vector) {
    Rcpp::NumericVector pv(p_mix_vec.begin(), p_mix_vec.end());
    resolved_priors["p_mix"] = pv;
  }
  if (m_nex_vector) {
    Rcpp::NumericVector nv(m_nex_vec.begin(), m_nex_vec.end());
    resolved_priors["m_nex"] = nv;
  }
  if (v_nex_vector) {
    Rcpp::NumericVector vv(v_nex_vec.begin(), v_nex_vec.end());
    resolved_priors["v_nex"] = vv;
  }

  Rcpp::List diagnostics = Rcpp::List::create(
    Rcpp::Named("n_obs") = n,
    Rcpp::Named("n_groups") = K,
    Rcpp::Named("n_covariates") = P,
    Rcpp::Named("n_censored") = arma::sum(1 - event),
    Rcpp::Named("n_observed") = arma::sum(event),
    Rcpp::Named("resolved_priors") = resolved_priors
  );

  return Rcpp::List::create(
    Rcpp::Named("draws") = draws_out,
    Rcpp::Named("priors") = priors,
    Rcpp::Named("iter") = iter,
    Rcpp::Named("warmup") = warmup,
    Rcpp::Named("chains") = chains,
    Rcpp::Named("diagnostics") = diagnostics
  );
}
