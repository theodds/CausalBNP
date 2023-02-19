#ifndef MATH_UTIS_H
#define MATH_UTILS_H

#include <RcppArmadillo.h>

double LogSumExp(arma::vec &x);
double log1m_inv_logit(double x);
double log_inv_logit(double x);

#endif
