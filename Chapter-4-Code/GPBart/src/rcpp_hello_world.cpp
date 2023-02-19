#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}
    
// [[Rcpp::export]]

arma::mat GetRankKernCat(const arma::mat& ranks, const arma::mat& cats) {
  double N = ranks.n_rows; 
  double P_qual = ranks.n_cols;
  double P_cat = cats.n_cols;
  double P = P_qual + P_cat;
  
  arma::mat out = arma::ones<arma::mat>(ranks.n_rows, ranks.n_rows);
  
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      for(int p = 0; p < P_qual; p++) {
        out(i,j) -= abs(ranks(i,p) - ranks(j,p)) / (N * P);
      }
      for(int p = 0; p < P_cat; p++) {
        out(i,j) -= (cats(i,p) == cats(j,p) ? 0.0 : 0.5 / P);
      }
    }
  }
  
  return out;
}

// [[Rcpp::export]]
arma::mat GetRankKern(const arma::mat& ranks) {
  double N = ranks.n_rows; 
  double P = ranks.n_cols;
  
  arma::mat out = arma::zeros<arma::mat>(ranks.n_rows, ranks.n_rows);
  
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      out(i,j) += 1;
      for(int p = 0; p < P; p++) {
        out(i,j) -= abs(ranks(i,p) - ranks(j,p)) / (N * P);
      }
    }
  }
  
  Rcout << "out[1,1] = " << out[1,1];
  return out;
}

// [[Rcpp::export]]
void test_const_ref( const arma::mat& m ){}

// [[Rcpp::export]]
double log_sum_exp(arma::vec x) {
  double m = x.max();
  return m + log(sum(exp(x-m)));
}
