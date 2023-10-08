#include <RcppArmadillo.h>
#include "MathUtils.h"
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec ParafacGcomp(int N_sim,
                  int J,
                  int j_0,
                  const arma::vec& omega,
                  const arma::vec& log_omega,
                  const arma::mat& gamma,
                  const arma::mat& log_gamma,
                  const arma::mat& log_1_m_gamma,
                  const arma::mat& beta,
                  const arma::mat& log_beta,
                  const arma::mat& log_1_m_beta,
                  const double& xi
                  ) {


  arma::vec Y_samp(N_sim);
  int K = omega.n_elem;

  int j_samp = j_0 - 1;

  for(int n = 0; n < N_sim; n++) {


    arma::uvec R(J);
    arma::uvec Y(J);
    int C = sample_class(omega);


    for(int j = 0; j < J; j++) {
      R(j) = R::rbinom(1.0, gamma(C,j));
      Y(j) = R::rbinom(1.0, beta(C,j));
    }

    if(R(j_samp) == 1) {
      // Y_samp(n) = Y(j_0);
      Y_samp(n) = beta(C,j_samp);
    } else {
      R(j_samp) = 1;
      Y(j_samp) = 1;
      arma::vec log_probs_success = log_omega;
      arma::vec log_probs_failure = log_omega;
      for(int k = 0; k < K; k++) {
        for(int j = 0; j < J; j++) {
          if(R(j) == 1) {
            log_probs_success(k) = log_probs_success(k) + log_gamma(k,j);
            if(Y(j) == 1) {
              log_probs_success(k) = log_probs_success(k) + log_beta(k,j);
            }
            else {
              log_probs_success(k) = log_probs_success(k) + log_1_m_beta(k,j);
            }
          }
          else {
            log_probs_success(k) = log_probs_success(k) + log_1_m_gamma(k,j);
          }
        }
        log_probs_failure(k) = log_probs_success(k) - log_beta(k,j_samp) + log_1_m_beta(k,j_samp);
      }
      arma::vec f(2);
      f(0) = LogSumExp(log_probs_success) + xi;
      f(1) = LogSumExp(log_probs_failure);
      double success_prob = exp(f(0) - LogSumExp(f));

      Y_samp(n) = success_prob;

    }

  }


  return Y_samp;
  // return List::create(Named("Y_samp") = Y_samp);

}
