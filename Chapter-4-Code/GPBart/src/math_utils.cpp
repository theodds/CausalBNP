#include "math_utils.h"

using namespace arma;
using namespace Rcpp;

double LogSumExp(vec &x){
  double m = x.max();
  int n = x.size();
  double out = 0;
  for(int i = 0; i < n; i++){
    out += exp(x(i) - m);
  }
  out = m + log(out);
  return(out);
}

double log1m_inv_logit(double x) {
  if(x > 0.0) 
    return -x - log1p(exp(-x));
  return -log1p(exp(x));
}

double log_inv_logit(double x) {
  if(x < 0.0)
    return x - log1p(exp(x));
  return -log1p(exp(-x));
}

