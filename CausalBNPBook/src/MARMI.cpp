#include <RcppArmadillo.h>
#include "MathUtils.h"
using namespace Rcpp;

// [[Rcpp::export]]
arma::umat MARMI(const arma::umat& Y,
                const arma::umat& R,
                const arma::vec& omega,
                const arma::vec& log_omega,
                const arma::mat& beta,
                const arma::mat& log_beta,
                const arma::mat& log_1_m_beta) {

  int N = Y.n_rows;
  int J = Y.n_cols;
  int K = log_omega.n_elem;
  arma::umat Y_out = Y;

  // Imputation for each observation
  for(int n = 0; n < N; n++) {

    // Compute class membership probability
    arma::vec logliks = log_omega;
    for(int k = 0; k < K; k++) {
      for(int j = 0; j < J; j++) {
        if(R(n,j) == 1) {
          logliks(k) = logliks(k) +
            (Y(n,j) == 1 ? log_beta(k,j) : log_1_m_beta(k,j));
        }
      }
    }

    // Sample latent class
    arma::vec class_probs = exp(logliks - LogSumExp(logliks));
    int C = sample_class(class_probs);

    // Impute the missing values
    for(int j = 0; j < J; j++) {
      if(R(n,j) == 0) {
        Y_out(n,j) = log(unif_rand()) < log_beta(C,j) ? 1 : 0;
      }
    }

  }

  return Y_out;

}

