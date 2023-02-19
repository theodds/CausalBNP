#include "logliks.h"


double sigma_inv_gamma_loglik::g(double x) {
  double tau = 1 / (x * x);
  return (
	  shape * log(rate) - R::lgammafn(shape) + (shape - 1) * log(tau) - rate * tau 
	  + 0.5 * N * log(tau) - 0.5 * tau * SSE
	  );
}

double sigma_log_normal_loglik::g(double x) {
  double tau = 1 / (x * x);
  double logx = log(x);
  double prior_weight = -logx - pow(logx - mu,2) / (2 * pow(sd, 2));
  double loglik_weight = 0.5 * N * log(tau) - 0.5 * tau * SSE;
  return(prior_weight + loglik_weight);
}

double sigma_log_normal_loglik::sample(double x) {
  return(uni_slice<sigma_log_normal_loglik>(x, *this, 1.0, INT_MAX, 0.0, 
					    ML_POSINF, ML_NAN));
}

double sigma_uniform_loglik::g(double x) {
  double tau = 1 / (x * x);
  return(0.5 * N * log(tau) - 0.5 * tau * SSE);
}

double sigma_uniform_loglik::sample(double x) {
  return(uni_slice<sigma_uniform_loglik>(x, *this, 1.0, INT_MAX, 0.0, 
					 upper, ML_NAN));
}

double mu_cauchy_normal::g(double mu) {
  double prior_weight = -log(1.0 + pow((mu - loc_mu) / scale_mu, 2));
  double loglik_weight = -0.5 * prec_Y * SSE_Y - 
    0.5 * N * prec_Y * pow(Y_bar - mu, 2);
  return(prior_weight + loglik_weight);
}

double mu_cauchy_normal::sample(double mu) {
  return(uni_slice<mu_cauchy_normal>(mu, *this, 1.0, INT_MAX, ML_NEGINF, 
				     ML_POSINF, ML_NAN));
}

double alpha_lognormal_dirichlet::g(double alpha) {
  double prior_weight = -log(alpha) 
    - 0.5 * prec_alpha * pow(log(alpha) - mu_alpha, 2);
  double loglik_weight = R::lgammafn(alpha) - R::lgammafn(alpha / K) + 
    alpha * log_V_bar;
  return(prior_weight + loglik_weight);
}

double alpha_lognormal_dirichlet::sample(double alpha) {
  return(uni_slice<alpha_lognormal_dirichlet>(alpha, *this, 1.0, INT_MAX, 
					      0.0, ML_POSINF, ML_NAN));
}

double intercept_single_pred_glm_cauchy_prior::g(double intercept) {
  double prior_weight = -log(1.0 + pow((intercept - loc_intercept) / scale_intercept, 2));
  double loglik_weight = 0;
  for(int n = 0; n < X.size(); n++) {
    if(Y(n) == 1) 
      loglik_weight += log_inv_logit(intercept + slope * X(n));
    else
      loglik_weight += log1m_inv_logit(intercept + slope * X(n));
  }
  return(prior_weight + loglik_weight);
}

double intercept_single_pred_glm_cauchy_prior::sample(double intercept) {
  return(uni_slice<intercept_single_pred_glm_cauchy_prior>(intercept, *this, 1.0, INT_MAX, 
							   ML_NEGINF, ML_POSINF, ML_NAN));
}

double slope_single_pred_glm_cauchy_prior::g(double slope) {
  double prior_weight = -log1p(pow((slope - loc_slope) / scale_slope, 2));
  double loglik_weight = 0.0;
  for(int n = 0; n < X.size(); n++) {
    if(Y(n) == 1)
      loglik_weight += log_inv_logit(intercept + slope * X(n));
    else
      loglik_weight += log1m_inv_logit(intercept + slope * X(n)); 
  }
  return(prior_weight + loglik_weight);
}

double slope_single_pred_glm_cauchy_prior::sample(double slope) {
  return(uni_slice<slope_single_pred_glm_cauchy_prior>(slope, *this, 1.0, INT_MAX, 
						       ML_NEGINF, ML_POSINF, ML_NAN));
}
