#include "Node.h"

int Node::NumBotNodes() {
  if(Bot) {
    return 1; 
  }
  else {
    return(Left->NumBotNodes() + Right->NumBotNodes());
  }
}

// arma::uvec Node::GetVec(arma::vec& x) {
//   // Should never run this unless something is wrong
//   if(Top == 0) return 9999 * arma::ones<arma::uvec>(3); 

//   arma::uvec out = arma::zeros<arma::uvec>(NumBot); 

//   int lab = LabelCov(arma::vec& x);
//   out(lab) = 1;

//   return(out); 
// }

int Node::LabelCov(arma::vec& x) {
  if(Bot) {
    return(Label);
  }
  else if(x(predictor) <= value) {
    return(Left->LabelCov(x));
  } else {
    return(Right->LabelCov(x));
  }
}

double LogSplitProb(const double& alpha, const double& beta, const int& depth) {
  return log(alpha) - beta * log(1 + depth); 
}

int DoSplit(const double& alpha, const double& beta, const int& depth) {
  int retval = 0; 
  double log_prob = LogSplitProb(alpha, beta, depth); 
  double log_cutoff = -exp_rand();

  if(log_cutoff <= log_prob) {
    retval = 1;
  }
  return retval;
}

int Sample(arma::vec& p) {
  int out = 0;
  double cs = 0.0; 
  double cutoff = unif_rand(); 
  for(int i = 0; i < p.size(); i++) {
    cs += p[i];
    if(cutoff <= cs) {
      out = i;
      break;
    }
  }
  return(out); 
}

void Node::InitTree(double& alpha, double& beta, arma::vec& cov_select) {
  GrowTree(alpha, beta, 0, this, cov_select); 
  LabelChild(0); 
}

void Node::GrowTree(double& alpha, double& beta, int depth, Node* parent, arma::vec& cov_select) {
  Depth = depth; 
  Top = (depth == 0) ? 1 : 0;
  Parent = parent; 

  int do_split = DoSplit(alpha, beta, depth); 
  if(do_split == 0) {
    Bot = 1; 
  } else {
    Left = new Node(); 
    Right = new Node(); 
    Right->GrowTree(alpha, beta, depth + 1, this, cov_select);
    Left->GrowTree(alpha, beta, depth + 1, this, cov_select);

    value = unif_rand(); 
    predictor = Sample(cov_select); 
    Bot = 0;
  }
}

Node::Node() {
  Top = 1;
  Bot = 1;
  Depth = 0;
}

int Node::LabelChild(int lab) {
  // Associates each leaf node to a label
  if(Bot) {
    Label = lab;
    return(lab + 1); 
  } else {
    lab = Left->LabelChild(lab); 
    lab = Right->LabelChild(lab); 
    return(lab); 
  }
}

// [[Rcpp::export]]
arma::mat TreeCompress(double& alpha, double& beta, arma::vec& cov_select, arma::mat& X) {
  Node* tree = new Node(); 
  tree->InitTree(alpha, beta, cov_select); 

  int num_bottom = tree->NumBotNodes(); 
  int N = X.n_rows;
  arma::mat X_star = arma::zeros<arma::mat>(N, num_bottom); 
  for(int i = 0; i < N; i++) {
    arma::vec x = arma::trans(X.row(i)); 
    int lab = tree->LabelCov(x); 
    X_star(i, lab) = 1.0;
  }

  return X_star;
}
