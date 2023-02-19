#ifndef LOGLIKS_H
#define LOGLIKS_H

#include <RcppArmadillo.h>
#include "math_utils.h"
#include "Slice.h"

// struct loglik {
//   Rcpp::NumericVector Y;
//   double nu;
//   double psi;
  
//   double g(double x, int param);
// };

/*THIS PRIOR RETURNS THINGS IN TERMS OF SIGMA, THE SD*/
/*ACTUALLY, ALL PRIORS ARE IN TERMS OF SIGMA*/

struct sigma_inv_gamma_loglik {
  double SSE;
  double shape;
  double rate;
  int N;
  
  double g(double x);
};

struct sigma_log_normal_loglik {
  
  double SSE;
  double mu;
  double sd;
  int N;

  double g(double x);
  double sample(double x);
};

struct sigma_uniform_loglik {
  double SSE;
  int N;
  double upper;

  double g(double x);
  double sample(double x);
};

/*If the data are heteroskeadstic, take prec_Y = 1, 
Y_bar = the weighted mean, SSE_Y the weighted SSE, 
and N = sum(precisions)*/
struct mu_cauchy_normal {
  double loc_mu;
  double scale_mu;
  double prec_Y;
  double SSE_Y;
  double Y_bar;
  double N;

  double g(double mu);
  double sample(double x);
};



struct alpha_lognormal_dirichlet {
  double mu_alpha;
  double prec_alpha;
  double log_V_bar;
  int K;

  double g(double alpha);
  double sample(double alpha);
};

struct intercept_single_pred_glm_cauchy_prior {
  arma::uvec Y;
  arma::vec X;
  double loc_intercept;
  double scale_intercept;
  double slope;

  double g(double intercept);
  double sample(double intercept);
};

struct slope_single_pred_glm_cauchy_prior {
  arma::uvec Y;
  arma::vec X;
  double loc_slope;
  double scale_slope;
  double intercept;

  double g(double slope);
  double sample(double slope);
};

#endif
