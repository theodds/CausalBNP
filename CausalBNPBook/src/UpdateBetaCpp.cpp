#include <RcppArmadillo.h>
#include "MathUtils.h"
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::mat UpdateBetaCpp(arma::mat success_counts, arma::mat failure_counts,
                        arma::vec col_shape_1, arma::vec col_shape_2) {

  int K = success_counts.n_rows;
  int J = success_counts.n_cols;

  arma::mat beta_out(K, J);

  for(int k = 0; k < K; k++) {
    for(int j = 0; j < J; j++) {
      beta_out(k,j) = R::rbeta(success_counts(k,j) + col_shape_1(j),
                               failure_counts(k,j) + col_shape_2(j));
    }
  }


  return beta_out;
}

// [[Rcpp::export]]
List UpdateSufficient(arma::umat Y, arma::umat R, arma::uvec C, int K) {

  int N = R.n_rows;
  int J = R.n_cols;

  arma::mat success_counts_Y = arma::zeros<arma::mat>(K, J);
  arma::mat success_counts_R = arma::zeros<arma::mat>(K, J);
  arma::mat failure_counts_Y = arma::zeros<arma::mat>(K, J);
  arma::mat failure_counts_R = arma::zeros<arma::mat>(K, J);
  arma::vec class_counts = arma::zeros<arma::vec>(K);

  for(int n = 0; n < N; n++) {
    int k = C(n) - 1;
    class_counts(k) = class_counts(k) + 1;
    for(int j = 0; j < J; j++) {
      if(R(n,j) == 1) {
        success_counts_R(k,j) = success_counts_R(k,j) + 1;
        if(Y(n,j) == 1) {
          success_counts_Y(k,j) = success_counts_Y(k,j) + 1;
        } else {
          failure_counts_Y(k,j) = failure_counts_Y(k,j) + 1;
        }
      } else {
        failure_counts_R(k,j) = failure_counts_R(k,j) + 1;
      }
    }
  }

  return List::create(Named("success_counts_Y") = success_counts_Y,
                      Named("success_counts_R") = success_counts_R,
                      Named("failure_counts_Y") = failure_counts_Y,
                      Named("failure_counts_R") = failure_counts_R,
                      Named("class_counts") = class_counts);

}






// [[Rcpp::export]]
List UpdateClassCpp(arma::umat Y,
                    arma::umat R,
                    arma::mat log_beta,
                    arma::mat log_1_m_beta,
                    arma::mat log_gamma,
                    arma::mat log_1_m_gamma,
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
          logliks(k) = logliks(k) + log_gamma(k,j);
          if(Y(n,j) == 1) {
            logliks(k) = logliks(k) + log_beta(k,j);
          }
          else {
            logliks(k) = logliks(k) + log_1_m_beta(k,j);
          }
        }
        else {
          logliks(k) = logliks(k) + log_1_m_gamma(k,j);
        }
      }
    }
    loglik_y(n) = LogSumExp(logliks);
    C(n) = sample_class_log(logliks);
  }

  return List::create(Named("C") = C, Named("loglik_data") = loglik_y);

}





