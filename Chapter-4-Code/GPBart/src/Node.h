#include <RcppArmadillo.h>
using namespace Rcpp;

class Rule {
 public: 

  int Var;
  double split; 
  int Right(const arma::vec& x); 
};

class Node {
 public: 
   
   Node(); 

  int Bot;
  int Top;
  int NumBot;
  int Label;
  int Depth;
  
  int predictor;
  double value;

  Node *Parent;
  Node *Left;
  Node *Right;

  int NumBotNodes();
  int LabelCov(arma::vec& x);
  void GrowTree(double& alpha, double& beta, int depth, Node* parent, arma::vec& cov_select);
  void InitTree(double& alpha, double& beta, arma::vec& cov_select);
  int LabelChild(int lab);
};

