#include "GPBart.h"

using namespace arma;

// [[Rcpp::export]]
Rcpp::List GPFit(const arma::vec& y_train,
                 const arma::mat& X_train,
                 const arma::mat& X_test,
                 const arma::mat& K_train,
                 const arma::mat& K_test,
                 const double sigma_sq,
                 const double sigma_sq_mu,
                 const double lambda,
                 const double rate_tau,
                 const double shape_tau,
                 const int nburn,
                 const int nthin,
                 const int nsave,
                 const int update_sigma_sq_mu,
                 const int update_lambda)
{
  // Initialize some stuff
  GPClass gp;
  gp.y           = y_train;
  gp.K           = K_train; 
  gp.sigma_sq    = sigma_sq;
  gp.tau         = 1.0 / sigma_sq;
  gp.sigma_sq_mu = sigma_sq_mu;
  gp.lambda      = lambda;
  gp.rate_tau    = rate_tau;
  gp.shape_tau   = shape_tau;

  // Some derived constants
  int N      = y_train.size(); 
  int P      = X_train.n_cols;
  int N_test = X_test.n_rows;

  // Make output
  mat f_test       = zeros<mat>(N_test, nsave); 
  mat vs           = zeros<mat>(N, nsave);
  vec lambdas      = zeros<vec>(nsave);
  vec sigma_sqs    = zeros<vec>(nsave);
  vec sigma_sq_mus = zeros<vec>(nsave); 
  vec evidences    = zeros<vec>(nsave); 

  // Run burn in
  vec v; 
  double evidence;
  double w_lambda = 1.0;
  double w_tau = sqrt(shape_tau) / rate_tau; 
  double w_sigma_mu = 1.0 / 8.0;
  for(int i = 0; i < nburn; i++) {
    Train(y_train, K_train, gp.sigma_sq, gp.sigma_sq_mu, gp.lambda, v, evidence); 
    if(update_lambda == 1)
      {
        gp.lambda      = uni_slice(gp.lambda, gp, w_lambda, 10000, 0.0, 1.0/0.0, gp.g(gp.lambda));
      }
    if(update_sigma_sq_mu == 1)
      {
        gp.sigma_sq_mu = uni_slice2(gp.sigma_sq_mu, gp, w_sigma_mu, 10000, 0.0, 1.0/0.0,
                                    gp.g2(gp.sigma_sq_mu));
      }
    gp.sigma_sq    = 1.0 / gp.tau;
    gp.sigma_sq_mu = uni_slice2(gp.sigma_sq_mu, gp, w_sigma_mu, 10000, 0.0, 1.0/0.0,
                                gp.g2(gp.sigma_sq_mu));
    Rcpp::Rcout << "Finishing Iteration" << i << std::endl;
  }

  for(int i = 0; i < nsave; i++) {
    for(int j = 0; j < nthin; j++) {
       
      if(update_lambda == 1)
        {
          gp.lambda      = uni_slice(gp.lambda, gp, w_lambda, 10000, 0.0, 1.0/0.0, gp.g(gp.lambda));
        }
      if(update_sigma_sq_mu == 1)
        {
          gp.sigma_sq_mu = uni_slice2(gp.sigma_sq_mu, gp, w_sigma_mu, 10000, 0.0, 1.0/0.0,
                                      gp.g2(gp.sigma_sq_mu));
        }
      gp.tau         = uni_slice1(gp.tau, gp, w_tau, 10000, 0.0, 1.0/0.0, gp.g1(gp.tau));
      gp.sigma_sq    = 1.0 / gp.tau;
      Train(y_train, K_train, gp.sigma_sq, gp.sigma_sq_mu, gp.lambda, v, evidence);
      Rcpp::Rcout << "Finishing Save" << i << std::endl;
    }
    lambdas(i)      = gp.lambda; 
    sigma_sqs(i)    = gp.sigma_sq;
    sigma_sq_mus(i) = gp.sigma_sq_mu; 
    evidences(i)    = evidence;
    f_test.col(i)   = gp.sigma_sq_mu * exp(-gp.lambda * K_test) * v;
    vs.col(i)       = v;
  }

  return Rcpp::List::create(
                            Rcpp::Named("lambda") = lambdas,
                            Rcpp::Named("sigma_sq") = sigma_sqs,
                            Rcpp::Named("sigma_sq_mu") = sigma_sq_mus,
                            Rcpp::Named("evidence") = evidences,
                            Rcpp::Named("f_test") = f_test,
                            Rcpp::Named("v") = vs
                            );

}


double GPClass::g(double x)
{
  return LogPriorLambda(x) + Evidence(y, K, sigma_sq, sigma_sq_mu, x);
}

double GPClass::g1(double x)
{
  return LogPriorTau(x, rate_tau, shape_tau) + Evidence(y, K, 1.0 / x, sigma_sq_mu, lambda);
}

double GPClass::g2(double x)
{
  return LogPriorSigmaMu(sqrt(x)) + Evidence(y, K, sigma_sq, x, lambda) - log(2) - 0.5 * log(x); 
}

double LogPriorTau(double x, double rate, double shape) {
  return LogDGamma(x, rate, shape);
}

double LogPriorLambda(double x) {
  return LogDExp(x, 1.0);
}

double LogPriorSigmaMu(double x) {
  return LogDExp(x, 8.0);
}

double LogDGamma(double x, double rate, double shape) {
  return rate * log(shape)  - lgamma(shape) + (rate - 1) * log(x) - rate * x;
}

double LogDExp(double x, double rate) {
  return log(rate) - rate * x;
}

double Evidence(const vec &y, const mat& K, const double sigma_sq,
                const double sigma_sq_mu,
                const double lambda)
{
  mat Sigma = sigma_sq_mu * exp(-lambda * K) + sigma_sq * eye(size(K));
  double log2pi = 1.8378770664093454835606594728112352797227949472755668;
  vec v = solve(Sigma, y);

  double evidence;
  double val, sign;
  log_det(val, sign, Sigma);

  evidence  = -0.5 * dot(y, v); 
  evidence -= 0.5 * y.size() * log2pi;
  evidence -= 0.5 * val;

  return evidence;
}

void Train(const vec &y, const mat& K, const double sigma_sq, const double sigma_sq_mu,
           const double lambda, vec& v, double& evidence)
{
  mat Sigma = sigma_sq_mu * exp(-lambda * K) + sigma_sq * eye(size(K));
  double log2pi = 1.8378770664093454835606594728112352797227949472755668;

  v = solve(Sigma, y);

  double val, sign;
  log_det(val, sign, Sigma);

  evidence = -0.5 * dot(y, v); 
  evidence -= 0.5 * y.size() * log2pi;
  evidence -= 0.5 * val;
}

vec Predict(const vec& v, const mat& K_test,
            const double sigma_sq_mu, const double lambda)
{
  return sigma_sq_mu * exp(-lambda * K_test) * v;
}
