#include <RcppArmadillo.h>
#include "MathUtils.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat NIPMI(const arma::mat& Y,
                const arma::umat& R,
                const arma::vec& omega,
                const arma::vec& log_omega,
                const arma::mat& gamma,
                const arma::mat& log_gamma,
                const arma::mat& log_1_m_gamma,
                const arma::mat& beta,
                const arma::mat& log_beta,
                const arma::mat& log_1_m_beta,
                const int j_0
                ) {

  int N = Y.n_rows;
  int J = Y.n_cols;
  int K = log_omega.n_elem;
  arma::mat Y_out = Y;
  int j_samp = j_0 - 1;

  // Imputation for each observation
  for(int n = 0; n < N; n++) {

    // Only need to do something if obs is missing
    if(R(n,j_samp) == 0) {

      // Compute class membership probability
      arma::vec logliks = log_omega;

      for(int k = 0; k < K; k++) {
        for(int j = 0; j < J; j++) {

          if(R(n,j) == 0) {
            logliks(k) = logliks(k) + log_1_m_gamma(k,j);
          }
          else {
            logliks(k) = logliks(k) + log_gamma(k,j); 
            if(Y(n,j) == 0) logliks(k) = logliks(k) + log_1_m_beta(k,j);
            if(Y(n,j) == 1) logliks(k) = logliks(k) + log_beta(k,j);
          }
        }

        logliks(k) = logliks(k) - log_1_m_gamma(k,j_samp) + log_gamma(k,j_samp);

      }

      // Sample latent class
      arma::vec class_probs = exp(logliks - LogSumExp(logliks));
      int C = sample_class(class_probs);

      // Impute the missing value
      Y_out(n,j_samp) = log(unif_rand()) < log_beta(C,j_samp) ? 1 : 0;
    }

  } 
  return Y_out;
}
