#include <RcppArmadillo.h>
#include "MathUtils.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat TLOMI(const arma::mat& Y,
                 const arma::umat& R,
                 const arma::vec& omega,
                 const arma::vec& log_omega,
                 const arma::mat& gamma,
                 const arma::mat& log_gamma,
                 const arma::mat& log_1_m_gamma,
                 const arma::mat& beta,
                 const arma::mat& log_beta,
                 const arma::mat& log_1_m_beta, 
                 const double xi,
                 const int j_0
) {
  
  int N = Y.n_rows;
  int J = Y.n_cols;
  int K = log_omega.n_elem;
  arma::mat Y_out = Y;
  int j_samp = j_0 - 1;
  
  // Imputation for each observation
  for(int n = 0; n < N; n++) {
    // Skip imputation if not needed
    if(R(n,j_samp) == 1) continue;
    
    // Compute class membership probability
    arma::vec logliks = log_omega;
    for(int k = 0; k < K; k++) {
      for(int j = 0; j < J; j++) {
        logliks(k) = logliks(k) + log_gamma(k,j);
        if(R(n,j) == 1) {
          logliks(k) = logliks(k) + (Y(n,j) == 1 ? log_beta(k,j) : log_1_m_beta(k,j));
        }
      }
    }
    
    // Compute untilted probability
    arma::vec class_probs = exp(logliks - LogSumExp(logliks));
    arma::vec success_probs = class_probs;
    for(int k = 0; k < K; k++) {
      success_probs(k) = success_probs(k) * beta(k,j_samp);
    }
    double success_prob = sum(success_probs);
    double fail_prob = 1 - success_prob;

    // Impute the missing value
    double tilted_success_prob = exp(xi) * success_prob / (exp(xi) * success_prob + fail_prob);
    Y_out(n,j_samp) = unif_rand() < tilted_success_prob ? 1 : 0;

  }
  
  return Y_out;
}
