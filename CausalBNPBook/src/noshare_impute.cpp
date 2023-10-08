#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

int sample_class(arma::vec &class_probabilities,
                 bool normalized = true, bool log_scale = false);
arma::mat rwish(arma::mat &Psi, double nu, bool doChol);
arma::vec mvrnorm(arma::vec &mu, arma::mat &Lambda, bool Precision = false);
arma::vec mvrnorm_chol(arma::vec &mu, arma::mat &R);
arma::vec rdirichlet(const arma::vec &alpha);
arma::mat mvrnorm(int N, arma::vec &mu, arma::mat &Sigma);

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


int sample_class(vec &class_probabilities,
		 bool normalized, bool log_scale) {
  double U = R::unif_rand();
  double foo = 0.0;
  if (normalized == true) {
    for(int k = 0; k < class_probabilities.size(); k++) {
      foo += class_probabilities(k);
      if(U < foo) {
	return(k);
      }
    }
  }
  if (normalized == false) {
    vec probs(class_probabilities.size());
    if(log_scale == true) {
      double logsumexp = LogSumExp(class_probabilities);
      for(int k = 0; k < class_probabilities.size(); k++) {
	probs(k) = exp(probs(k) - logsumexp);
      }
    }
    else {
      probs = class_probabilities / sum(class_probabilities);
    }
    for(int k = 0; k < class_probabilities.size(); k++) {
      foo += probs(k);
      if(U < foo) {
	return(k);
      }
    }
  }
}

mat rwish(mat &Psi, double nu, bool doChol) {
  mat S = Psi;
  if(doChol) {
    try {
      S = chol(S);
    }
    catch(...) {
      // Rcout << "Failed in Rwish doing Cholesky" << std::endl; //Debug
      // printMat(S);
      throw;
    }
  }
  int J = Psi.n_rows;
  mat B; B.zeros(J, J);
  for(int i = 0; i < J; i++){
    double df = nu - i;
    vec X(df); X.randn();
    B(i, i) = sqrt(as_scalar(trans(X) * X));
    if(i > 0) {
      vec foo(i); foo.randn();
      for(int j = 0; j < i; j++){
	B(i, j) = foo(j);
      }
    }
  }
  return(trans(S) * B * trans(B) * S);
}

vec mvrnorm(vec &mu, mat &Lambda, bool Precision) {
  vec X; X.randn(mu.size());
  mat Sigma;
  if(Precision == true) {Sigma = inv((Lambda));} else {Sigma = Lambda;}
  mat R = chol(Sigma);
  return(trans(R) * X + mu);
}
vec mvrnorm_chol(vec &mu, mat &R) {
  vec X; X.randn(mu.size());
  return(trans(R) * X + mu);
}

mat mvrnorm(int N, vec &mu, mat &Sigma) {
  int J = mu.size();
  mat out(N, J); out.zeros(N, J);
  mat R = chol(Sigma);
  for(int n = 0; n < N; n++) {
    out.row(n) = trans(mvrnorm_chol(mu, R));
  }
  return(out);
}

vec rdirichlet(const vec &alpha) {
  vec result; result.zeros(alpha.size());
  for(int j = 0; j < alpha.size(); j++) {
    result(j) = R::rgamma(alpha(j), 1);

  }
  result = result / sum(result);
  return(result);
}

RcppExport SEXP noshare_impute(SEXP sa, SEXP ssigma, SEXP sphi,
			     SEXP szeta, SEXP sgamma, SEXP sxi, SEXP snum_sim,
			     SEXP snum_iter, SEXP sJ, SEXP sK,
			     SEXP ssens_param) {

  RNGScope scope;

  NumericVector ra(sa);
  NumericVector rsigma(ssigma);
  NumericVector rphi(sphi);
  NumericVector rzeta(szeta);
  NumericVector rgamma(sgamma);
  NumericVector rxi(sxi);
  NumericMatrix rsens_param(ssens_param);
  int num_sim = as<int>(snum_sim);
  int num_iter = as<int>(snum_iter);
  int J = as<int>(sJ);
  int K = as<int>(sK);
  mat sens_param(rsens_param.begin(), num_iter, J);

  vec change_from_baseline = zeros<vec>(num_iter);
  vec change_from_baseline_var = zeros<vec>(num_iter);
  mat out; out.zeros(num_iter, J);
  mat out_var; out_var.zeros(num_iter, J);
  mat out_drop_probs; out_drop_probs.zeros(num_iter, J);
  mat out_haz; out_haz.zeros(num_iter, J);
  mat out_obs_mean; out_obs_mean.zeros(num_iter, J);
  mat out_cum_drop_probs; out_cum_drop_probs.zeros(num_iter, J);
  for(int i = 0; i < num_iter; i++) {
    mat Y; Y.zeros(num_sim, J);
    uvec S_sim; S_sim.zeros(num_sim);
    mat a(ra.begin() + K * J * i, K, J, false);
    mat sigma(rsigma.begin() + K * J * i, K, J, false);
    cube phi(rphi.begin() + K * J * J * i, K, J, J, false);
    mat zeta(rzeta.begin() + K * J * i, K, J, false);
    mat gamma(rgamma.begin() + K * J * i, K, J, false);
    vec xi(rxi.begin() + K * i, K, false);
    for(int n = 0; n < num_sim; n++) {

      // Impute Initial Sequence
      int C = sample_class(xi);
      int S = 0;
      double m = a(C,0);
      Y(n,0) = R::rnorm(m, sigma(C,0));
      double p = exp(log1m_inv_logit(zeta(C,0) + gamma(C,0) * Y(n,0)));
      int proceed = R::rbinom(1, p);

      int j = 1;
      while(proceed == 1 && j < J) {
      	m = a(C,j);
      	for(int l = 0; l < j; l++) m += phi(C,j,l) * Y(n,l);
      	Y(n,j) = R::rnorm(m, sigma(C, j));
      	p = exp(log1m_inv_logit(zeta(C,j) + gamma(C,j) * Y(n,j)));
      	proceed = R::rbinom(1,p);
      	j++;
      }
      S = j;
      S_sim(n) = S - 1;
      // Finish off the sequence
      // Step 1: get cluster logliks for ge j and == j-1
      // Step 2: Impute Y according to loglik_ge_j
      // Step 3: If j immediately follows the observed piece,
      //         transform, otherwise
      //         transform only with a specified probability.
      while(j < J) {
      	vec loglik_ge_j; loglik_ge_j.zeros(K);
      	vec loglik_jm1; loglik_jm1.zeros(K);
      	for(int k = 0; k < K; k++) {
      	  double m = a(k,0);
      	  loglik_ge_j(k) = log(xi(k));
      	  loglik_ge_j(k) += R::dnorm(Y(n,0), m, sigma(k,0), 1);
      	  // P(ge 0 + 1)
      	  loglik_ge_j(k) += log1m_inv_logit(zeta(k,0) + gamma(k,0) * Y(n,0));
      	  for(int l = 1; l < j; l++) {
      	    m = a(k,l);
      	    for(int w = 0; w < l; w++) m += phi(k,l,w) * Y(n,w);
      	    loglik_ge_j(k) += R::dnorm(Y(n,l), m, sigma(k,l), 1);
      	    // P(ge l + 1)
      	    loglik_ge_j(k) += log1m_inv_logit(zeta(k,l) + gamma(k,l) * Y(n,l));
      	  }
      	  loglik_jm1(k) = loglik_ge_j(k) -
      	    log1m_inv_logit(zeta(k,j-1) +  gamma(k,j-1)  * Y(n,j-1)) +
      	    log_inv_logit(zeta(k,j-1) + gamma(k,j-1) * Y(n, j-1));
      	}
      	double p_ge_j = LogSumExp(loglik_ge_j);
      	double p_jm1 = LogSumExp(loglik_jm1);
      	vec probs = exp(loglik_ge_j - p_ge_j);
      	int C = sample_class(probs);
      	double m = a(C,j);
      	for(int w = 0; w < j; w++) m += phi(C,j,w) * Y(n,w);
      	Y(n,j) = R::rnorm(m, sigma(C,j));
      	double trans_prob = 1.0 / (1.0 + exp(p_ge_j - p_jm1));
      	if(j == S) trans_prob = 1;
      	if(R::runif(0,1) < trans_prob) Y(n,j) += sens_param(i, j);
      	j++;
      }
    }
    out.row(i) = mean(Y);
    out_var.row(i) = var(Y) / num_sim;
    change_from_baseline(i) = mean(Y.col(J-1) - Y.col(0));
    change_from_baseline_var(i) = var(Y.col(J-1) - Y.col(0)) / num_sim;
    for(int j = 0; j < J; j++) {
      uvec idx_obs = find(S_sim >= j);
      vec tmp_Y = Y.col(j);
      out_obs_mean(i,j) = mean(tmp_Y.elem(idx_obs));
      out_drop_probs(i,j) = (double)(sum(S_sim == j)) / (double)(num_sim);
      if(j > 0) {
	     out_haz(i,j) = (double)(sum(S_sim == j - 1)) /
	       (double)(sum(S_sim >= j - 1));
      }
      out_cum_drop_probs(i,j) += out_drop_probs(i,j);
      if(j > 0)
	     out_cum_drop_probs(i,j) += out_cum_drop_probs(i, j-1);
    }
    if(i % 100 == 0)
      Rcout << "\rFinishing Iteration " << i << "        ";
  }
  Rcout << std::endl;
  List true_out;
  true_out["means"] = out;
  true_out["var_of_mean"] = out_var;
  true_out["obs.means"] = out_obs_mean;
  true_out["drop.probs"] = out_drop_probs;
  true_out["haz"] = out_haz;
  true_out["cum.drop.probs"] = out_cum_drop_probs;
  true_out["change_from_baseline"] = change_from_baseline;
  true_out["change_from_baseline_var"] = change_from_baseline_var;
  return(wrap(true_out));

}
