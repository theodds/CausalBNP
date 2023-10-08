#include <RcppArmadillo.h>

double LogSumExp(arma::vec x);
int sample_class_log(const arma::vec &log_probs);
int sample_class(const arma::vec &probs);
