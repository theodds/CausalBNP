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
double LogSumExp(arma::vec x) {
  return x.max() + log(sum(exp(x - x.max())));
}

int sample_class(const arma::vec &probs) {
  double U = R::unif_rand();
  double foo = 0.0;
  int K = probs.size();

  // Sample
  for(int k = 0; k < K; k++) {
    foo += probs(k);
    if(U < foo) {
      return(k);
    }
  }

  return K - 1;

}

int sample_class_log(const arma::vec &log_probs) {
  double U = R::unif_rand();
  double foo = 0.0;
  int K = log_probs.size();
  arma::vec probs(K);

  // Normalize
  double logsumexp = LogSumExp(log_probs);
  for(int k = 0; k < K; k++) {
    probs(k) = exp(log_probs(k) - logsumexp);
  }

  // Sample
  for(int k = 0; k < K; k++) {
    foo += probs(k);
    if(U < foo) {
      return(k);
    }
  }

  return K - 1;

}
