#ifndef GPBART_H
#define GPBART_H

#include <RcppArmadillo.h>
#include "Slice.h"

struct GPClass
{
  arma::vec y; 
  arma::mat K; 
  double sigma_sq;
  double sigma_sq_mu;
  double lambda;
  double tau; 
  double rate_tau;
  double shape_tau;

  double g(double lambda);
  double g1(double sigma_sq);
  double g2(double sigma_sq_mu);

  void UpdateLambda();
  void UpdateSigmaSq();
  void UpdateSigmaSqMu();
};

double LogPriorTau(double x, double rate, double shape);

double LogPriorLambda(double x); 

double LogPriorSigmaMu(double x);

double LogDGamma(double x, double rate, double shape);

double LogDExp(double x, double rate);

double Evidence(const arma::vec &y, const arma::mat& K, const double sigma_sq,
                const double sigma_sq_mu,
                const double lambda);

void Train(const arma::vec &y, const arma::mat& K, const double sigma_sq, const double sigma_sq_mu,
           const double lambda, arma::vec& v, double& evidence);

arma::vec Predict(const arma::vec& v, const arma::mat& K_test,
            const double sigma_sq_mu, const double lambda);

#endif
