#include <RcppArmadillo.h>
#include "MathUtils.h"
using namespace Rcpp;

// [[Rcpp::export]]
List UpdateClassMARCpp(arma::umat Y,
                       arma::umat R,
                       arma::mat log_beta,
                       arma::mat log_1_m_beta,
                       arma::vec log_omega) {

  int N = Y.n_rows;
  int J = Y.n_cols;
  int K = log_omega.n_elem;

  arma::uvec C = arma::zeros<arma::uvec>(N);
  arma::vec loglik_y(N);

  // Get matrix of log likelihoods
  for(int n = 0; n < N; n++) {
    arma::vec logliks = arma::zeros<arma::vec>(K);
    for(int k = 0; k < K; k++) {
      logliks(k) = log_omega(k);
      for(int j = 0; j < J; j++) {
        if(R(n,j) == 1) {
          if(Y(n,j) == 1) {
            logliks(k) = logliks(k) + log_beta(k,j);
          }
          else {
            logliks(k) = logliks(k) + log_1_m_beta(k,j);
          }
        }
      }
    }
    loglik_y(n) = LogSumExp(logliks);
    C(n) = sample_class_log(logliks);
  }

  return List::create(Named("C") = C, Named("loglik_data") = loglik_y);

}

