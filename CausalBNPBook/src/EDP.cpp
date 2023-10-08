#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <Rmath.h>
//#include "utilFunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

double expit(double x) {
  return( ( 1 / (1 + exp ( - x ) ) ) );
}


int rmultinomF(vec const& p) {
  vec csp = cumsum(p/sum(p));
  double rnd = runif(1)[0];
  int res = 0;
  int psize = p.size();

  for(int i = 0; i < psize; i++) {
    if(rnd>csp(i)) res = res+1;
  }

  return(res+1);
}

// taken from: http://gallery.rcpp.org/articles/simulate-multivariate-normal/
vec mvrnorm(vec mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(1, ncols);
  return mu + trans(Y * chol(sigma));
}

double rinvchisq(double df, double scale) {
  double res = as<double>(rchisq(1,df));
  return( (df*scale)/res );
}

//adapted from dmvnorm function from mvtnorm package in R
double dmvn(vec x, vec mu, mat sig, int logt) {
  mat dec = chol(sig);
  vec tmp = solve(dec,x-mu); //check
  double rss = sum(tmp%tmp);

  double logretval = -sum( log (dec.diag() ) ) - 0.5 * sig.n_cols * log (2 * M_PI) - 0.5 * rss;
  if(logt==0)
    logretval = exp(logretval);
  return logretval;
}

double updatevar(double x, double nu0, double tau0, double c0, double mu0) { //how to declare vars like nu0 and tau0 and c0??
  double varx;
  //int xsize = x.size();
  int xsize=1;
  double newdf = nu0 + xsize;
  if ( xsize==1 ) { varx = 0; }
  else if ( xsize > 1 ) { varx = var(x); }
  double numer = nu0*tau0+(xsize-1)*varx+(c0*xsize/(c0+xsize))*pow((x-mu0),2); //change x to mean(x) if using vector for x
  double newval = rinvchisq(newdf,numer/newdf);
  return newval;
}

double updatemean(double x, double tau, double c0, double mu0) {
  //int xsize = x.size();
  int xsize = 1;
  double newvar = 1 / (c0/tau+xsize/tau);
  double newmean = (mu0 * c0 / tau + x * xsize / tau) * newvar; //change x to mean(x) if using vector
  double newval = rnorm(1,newmean,sqrt(newvar))[0];
  return newval;
}

vec newbetafunction(vec curbet, vec x, int y, vec betainit, mat diagbetacov0) {
  int dimbet = curbet.size();
  mat propvar(dimbet,dimbet,fill::eye);
  vec proposed = mvrnorm(curbet, 0.01*propvar);

  //careful this is a sum in original
  double loglikenew = R::dbinom(y, 1, expit(dot(x,proposed)), 1) +
    dmvn(proposed,betainit,diagbetacov0,1) +
    dmvn(curbet,proposed,0.01*propvar,1);

  double loglikeold = R::dbinom(y, 1, expit(dot(x,curbet)), 1) +
    dmvn(curbet,betainit,diagbetacov0,1) +
    dmvn(proposed,curbet,0.01*propvar,1);

  double ans = std::min(exp(loglikenew-loglikeold),1.0);

  if(runif(1)[0] < ans)
    curbet = proposed;
  return curbet;
}

double newregvar(vec x, double y, double betaa0, double betab0, vec beta0) {
  double allp = x.size();
  double an = betaa0 + 0.5;
  mat prec0(allp,allp,fill::eye);
  mat newprec =  x*x.t() + prec0;
  vec betan = newprec.i() * (prec0*beta0 + x*y);
  mat part1 = beta0.t()*prec0*beta0;
  mat part2 = betan.t()*newprec*betan;
  if(part1.n_rows > 1 || part1.n_cols > 1 || part2.n_rows > 1 || part2.n_cols > 1) {
    Rprintf("Something wrong in newregvar function");
  }
  double bn = betab0 + 0.5*(y*y + part1(0,0) + part2(0,0));
  double outvar = 1/(R::rgamma(an,bn));
  return outvar;
}

vec newbet(vec x, double y, double sig, vec beta0) {
  double allp = x.size();
  mat prec0(allp,allp,fill::eye);
  mat newprec = x* x.t() + prec0;
  vec betan = newprec.i() * (prec0*beta0 + x * y);
  vec bet = mvrnorm(betan,sig*newprec.i());
  return bet;
}




// [[Rcpp::export]]
List cluster(arma::vec y,
             arma::mat X,
             arma::mat matX,
             arma::ivec Sy,
             arma::ivec Sx,
             arma::mat betaY,
             arma::mat xPiPars,
             arma::mat xMuPars,
             arma::mat xSigPars,
             double alphapsi,
             double alphatheta,
             arma::vec h0y,
             arma::vec h0i,
             arma::imat uniqueS,
             double c0,
             double mu0,
             double nu0,
             double tau,
             double a0,
             double b0,
             arma::vec betainit,
             arma::mat diagbetacov0,
             int p1,
             int ptx,
             int p2) {
  /*
   * y is nx1 vector of outcomes
   * x is nxp matrix of covariates
   * Sy is nx1 vector of Y cluster memberships
   * Sx is nx1 vector of X cluster memberships
   * betaY is kxp matrix of coef (k is number of Y clusters)
   * xPiPars is k*xp1 matrix of coef (k* is total # of clusters, p1 is # binary)
   * xMuPars is k*xp2 matrix of coef (p2 is # continuous covs)
   * xsigPars is k*xp2 matrix of coef
   * alphapsi is the current value of alpha_{psi}
   * alphatheta is the current value of alpha_{theta}
   * h0y is nx1 vector
   * h0i is nx1 vector
   * uniqueS is k*x2 matrix of ordered unique observations of Sy and Sx
   * indS is an nx1 indicator of which row of uniqueS everyone is in
   */

  int n = X.n_rows, p = X.n_cols;
  uvec indY, indX; // to contain vector of obs for cluster memberships

  // how many unique clusters
  int numY; //# unique Y clusters
  ivec uniqueY; //indicators for unique Y clusters (sorted)
  ivec numX; // one element for each Y cluster
  int numTotalCluster;

  // initialize some dummy vars to store things
  uvec ind_dummy, ind_dummy2;
  int num_dummy, num_dummy2;
  int count;
  double likeregy, prodx, prodx2;

  // containers for output
  int newCluster;

  rowvec newpipars(p1+ptx);
  rowvec newmupars(p2);
  rowvec newsigpars(p2);

  vec betadraw;
  vec newbeta;


  // LOOP THROUGH EVERY PERSON AND CHANGE CLUSTER MEMBERSHIPS
  for(int i = 0; i < n; i++) {

    //Rprintf("Loop: %d\n",i);
    // check if ith person is the lone person in his cluster
    ind_dummy = find(Sy == Sy(i) && Sx == Sx(i));
    num_dummy = ind_dummy.size();
    //Rprintf("Cluster Y,X,#:%d,%d,%d\n",Sy(i),Sx(i),num_dummy);

    if(num_dummy==1) { //if lone person in X-Y cluster
      //DELETE ASSOCIATED COEFFICIENTS IN Y AND X CLUSTER

      ind_dummy2 = find(Sy == Sy(i)); //check if only person in Y cluster too
      num_dummy2 = ind_dummy2.size();

      //Rprintf("num_dummy2: %d\n",num_dummy2);
      //Rprintf("Sy(i): %d\n", Sy(i));
      //Rcout << Sy.t() << std::endl;
      //Rcout << ind_dummy2.t() << std::endl;
      //Rcout << ind_dummy2.size() << std::endl;
      // delete Y coef if only one in Y cluster
      if(num_dummy2==1)
        betaY.shed_row(Sy(i)-1);

      // delete X coef
      //should find row in uniqueS that corresponds to person i
      ind_dummy = find(uniqueS.col(0)==Sy(i) && uniqueS.col(1)==Sx(i));

      xPiPars.shed_row(ind_dummy(0)); // CHECK --- NEEDED TO DO THIS BECAUSE ind_dummy is uvec
      xMuPars.shed_row(ind_dummy(0));
      xSigPars.shed_row(ind_dummy(0));

      //relabel X cluster
      for(int j = 0; j < Sx.size(); j++) {
        if(Sy(j) == Sy(i) && Sx(j) > Sx(i))
          Sx(j) = Sx(j) - 1;
      }

      for(int j = 0; j < uniqueS.n_rows; j++) {
        if(uniqueS(j,0) == Sy(i) && uniqueS(j,1) > Sx(i))
          uniqueS(j,1) = uniqueS(j,1) - 1;
      }

      //relabel Y cluster (if needed)
      if(num_dummy2==1) {
        for(int j = 0; j < Sy.size(); j++) {
          if(Sy(j) > Sy(i))
            Sy(j) = Sy(j) - 1;
        }

        for(int j = 0; j < uniqueS.n_rows; j++) {
          if(uniqueS(j,0) > Sy(i))
            uniqueS(j,0) = uniqueS(j,0) - 1;
        }
      }

      uniqueS.shed_row(ind_dummy(0)); //get rid of row
    }

    // NEED TO DELETE ROW OF Sy and Sx??? Think so
    Sy.shed_row(i); Sx.shed_row(i);


    // recalculate number of unique clusters
    numY = betaY.n_rows;
    numTotalCluster = xMuPars.n_rows;
    int totalposs = numY+numTotalCluster+1;
    //Rprintf("Y clusters, total clusters, Total possible clusters: %d,%d,%d\n",numY,numTotalCluster,totalposs);
    vec probs(totalposs);
    count=0;

    int njwoi, nljwoi; //counts for # in appropriate Y and X cluster,
    //excluding the ith person

    int numXj; //# of X clusters within each Y cluster


    for(int j = 0; j < numY; j++) {
      //FILL IN PROBS FOR EXISTING CLUSTERS

      //get count of number of X clusters within jth Y cluster
      ind_dummy = find(uniqueS.col(0) == (j+1));
      numXj = ind_dummy.size();
      //Rprintf("NumX cluster within Y cluster: %d\n",numXj);

      //get number of subjects within jth cluster
      ind_dummy = find(Sy == (j+1));
      njwoi = ind_dummy.size();

      //Rprintf("njwoi: %d\n",njwoi);
      // likelihood for each existing Y cluster
      likeregy = R::dbinom(y(i), 1 , expit( dot( matX.row(i),betaY.row(j) ) ) , 0);

      //Rprintf("likelihood for Y: %.2f\n",likeregy);
      for(int k = 0; k < numXj; k++) {
        prodx = 1;
        prodx2 = 1;

        ind_dummy = find(Sy == (j+1) && Sx == (k+1));
        nljwoi = ind_dummy.size();
        //Rprintf("nljwoi: %d\n",nljwoi);

        //likelihood for binary covariates
        for(int l = 0; l < ptx+p1; l++) {
          prodx = prodx*R::dbinom(X(i,l), 1, xPiPars(count, l), 0);
        }

        //Rprintf("Prodx: %.2f\n",prodx);

        //likelihood for continuous covariates
        for(int l = 0; l < p2; l++) {
          prodx2 = prodx2*R::dnorm(X(i,ptx+p1+l), xMuPars(count,l), sqrt(xSigPars(count,l)), 0 );
        }

        //Rprintf("Prodx2: %.2f\n",prodx2);

        probs(count) = ((njwoi*nljwoi)/(njwoi+alphapsi))*likeregy*prodx*prodx2;
        //Rprintf("\t current count and prob: %d, %.2f\n",count,probs(count));
        count++;
      }
    }


    for(int j = 0; j < numY; j++) {
      //FILL IN PROBS FOR NEW X CLUSTERS IN EXISTING Y CLUSTERS

      ind_dummy = find(Sy == (j+1));
      njwoi = ind_dummy.size();

      // likelihood for each existing Y cluster
      likeregy = R::dbinom(y(i), 1 , expit( dot( matX.row(i),betaY.row(j) ) ) , 0);
      probs(numTotalCluster + j) = ( (njwoi*alphapsi) / (njwoi+alphapsi) ) * likeregy*h0i(i);
      //Rprintf("New X cluster\n");
      //Rprintf("\t current count and prob: %d, %.2f\n",numTotalCluster+j,probs(numTotalCluster+j));
    }

    probs(numY+numTotalCluster) = alphatheta*h0y(i)*h0i(i); //prob for new Y cluster
    //Rprintf("New Y cluster prob\n");
    //Rprintf("\t current count and prob: %d, %.2f\n",numY+numTotalCluster,probs(numY+numTotalCluster));

    //USE MULTINOMIAL DISTRIBUTION TO CHOOSE NEW CLUSTER
    newCluster = rmultinomF(probs);
    probs.zeros();

    //Rprintf("The new cluster is: %d\n", newCluster);
    //need to map this integer to one of the clusters (or a new cluster)

    if(newCluster<=numTotalCluster)
    {
      Sy.insert_rows(i,1);
      Sy(i) = uniqueS(newCluster-1,0);
      Sx.insert_rows(i,1);
      Sx(i) = uniqueS(newCluster-1,1);
    }
    else
    {
      //find out whether this is a new Y cluster or X cluster
      if(newCluster == (numTotalCluster+numY+1))
      {

        Sx.insert_rows(i,1);
        Sx(i) = 1;

        // need functions rmvn, newbetafunction, updatevar, updatemean

        betadraw = mvrnorm(betainit, diagbetacov0);
        newbeta = newbetafunction(betadraw,trans(matX.row(i)),y(i),betainit,diagbetacov0);
        betaY.insert_rows(numY,newbeta.t());

        for(int j = 0; j<(p1+ptx); j++) {
          newpipars(j) = R::rbeta(X(i,j)+a0,1-X(i,j)+b0);
        }

        for(int j = 0; j < p2; j++) {
          newsigpars(j) = updatevar(X(i,(p1+ptx+j)),nu0,tau,c0,mu0);
          newmupars(j) = updatemean(X(i,(p1+ptx+j)),tau,c0,mu0);
        }
        xPiPars.insert_rows(numTotalCluster,newpipars);
        xSigPars.insert_rows(numTotalCluster,newsigpars);
        xMuPars.insert_rows(numTotalCluster,newmupars);


        Sy.insert_rows(i,1);
        Sy(i) = Sy.max()+1;
        uniqueS.insert_rows(numTotalCluster,1);
        uniqueS(numTotalCluster,0) = Sy(i);
        uniqueS(numTotalCluster,1) = Sx(i);
      }
      else
      {
        //Y cluster should be numCluster - numTotalCluster
        //X cluster should be max +1 of X clusters within Y

        Sy.insert_rows(i,1);
        Sy(i) = newCluster - numTotalCluster;
        ind_dummy = find(uniqueS.col(0) == Sy(i));
        Sx.insert_rows(i,1);
        Sx(i) = ind_dummy.size() + 1;

        // update parameters
        for(int j = 0; j<(p1+ptx); j++) {
          newpipars(j) = R::rbeta(X(i,j)+a0,1-X(i,j)+b0);
        }

        for(int j = 0; j < p2; j++) {
          newsigpars(j) = updatevar(X(i,(p1+ptx+j)),nu0,tau,c0,mu0);
          newmupars(j) = updatemean(X(i,(p1+ptx+j)),tau,c0,mu0);
        }

        ind_dummy = find(uniqueS.col(0) <= Sy(i));
        num_dummy = ind_dummy.size();

        xPiPars.insert_rows(num_dummy,newpipars);
        xSigPars.insert_rows(num_dummy,newsigpars);
        xMuPars.insert_rows(num_dummy,newmupars);

        uniqueS.insert_rows(num_dummy,1);
        uniqueS(num_dummy,0) = Sy(i);
        uniqueS(num_dummy,1) = Sx(i);
      }
    }


    // as necessary add back into Sy, Sx, uniqueS, etc...
    //Rprintf("New cluster for person i: %d, %d\n\n",Sy(i),Sx(i));
    //Rcout << uniqueS << std::endl;
    //Rcout << betaY << std::endl;
    //Rprintf("Max Sy, Max uniqueS: %d,%d\n",Sy.max(),uniqueS.col(0).max());
  }



  return List::create(_["Sy"]  = Sy,
                      _["Sx"] = Sx,
                      _["uniqueS"] = uniqueS,
                      _["beta"]    = betaY,
                      _["pipars"]  = xPiPars,
                      _["mupars"]  = xMuPars,
                      _["sigpars"] = xSigPars);
}

// [[Rcpp::export]]
List cluster_continuous(arma::vec y,
                        arma::mat X,
                        arma::mat matX,
                        arma::ivec Sy,
                        arma::ivec Sx,
                        arma::mat betaY,
                        arma::vec sig2,
                        arma::mat xPiPars,
                        arma::mat xMuPars,
                        arma::mat xSigPars,
                        double alphapsi,
                        double alphatheta,
                        arma::vec h0y,
                        arma::vec h0i,
                        arma::imat uniqueS,
	                      double c0,
	                      double mu0,
	                      double nu0,
	                      double tau,
	                      double a0,
	                      double b0,
	                      double betaa0,
	                      double betab0,
	                      arma::vec betainit,
	                      arma::mat diagbetacov0,
	                      int p1,
	                      int ptx,
	                      int p2) {
  /*
   * y is nx1 vector of outcomes
   * x is nxp matrix of covariates
   * Sy is nx1 vector of Y cluster memberships
   * Sx is nx1 vector of X cluster memberships
   * betaY is kxp matrix of coef (k is number of Y clusters)
   * xPiPars is k*xp1 matrix of coef (k* is total # of clusters, p1 is # binary)
   * xMuPars is k*xp2 matrix of coef (p2 is # continuous covs)
   * xsigPars is k*xp2 matrix of coef
   * alphapsi is the current value of alpha_{psi}
   * alphatheta is the current value of alpha_{theta}
   * h0y is nx1 vector
   * h0i is nx1 vector
   * uniqueS is k*x2 matrix of ordered unique observations of Sy and Sx
   * indS is an nx1 indicator of which row of uniqueS everyone is in
   */

  int n = X.n_rows, p = X.n_cols;
  uvec indY, indX; // to contain vector of obs for cluster memberships

  // how many unique clusters
  int numY; //# unique Y clusters
  ivec uniqueY; //indicators for unique Y clusters (sorted)
  ivec numX; // one element for each Y cluster
  int numTotalCluster;

  // initialize some dummy vars to store things
  uvec ind_dummy, ind_dummy2;
  int num_dummy, num_dummy2;
  int count;
 double likeregy, prodx, prodx2;

  // containers for output
  int newCluster;

  rowvec newpipars(p1+ptx);
  rowvec newmupars(p2);
  rowvec newsigpars(p2);

  vec betadraw;
  vec newbeta;
  double newsig;


  // LOOP THROUGH EVERY PERSON AND CHANGE CLUSTER MEMBERSHIPS
  for(int i = 0; i < n; i++) {

    //Rprintf("Loop: %d\n",i);
    // check if ith person is the lone person in his cluster
    ind_dummy = find(Sy == Sy(i) && Sx == Sx(i));
    num_dummy = ind_dummy.size();
    //Rprintf("Cluster Y,X,#:%d,%d,%d\n",Sy(i),Sx(i),num_dummy);

    if(num_dummy==1) { //if lone person in X-Y cluster
      //DELETE ASSOCIATED COEFFICIENTS IN Y AND X CLUSTER

      ind_dummy2 = find(Sy == Sy(i)); //check if only person in Y cluster too
      num_dummy2 = ind_dummy2.size();

      //Rprintf("num_dummy2: %d\n",num_dummy2);
      //Rprintf("Sy(i): %d\n", Sy(i));
      //Rcout << Sy.t() << std::endl;
      //Rcout << ind_dummy2.t() << std::endl;
      //Rcout << ind_dummy2.size() << std::endl;
      // delete Y coef if only one in Y cluster
      if(num_dummy2==1) {
	betaY.shed_row(Sy(i)-1);
	sig2.shed_row(Sy(i)-1);  // NEW FOR CONTINUOUS
      }
      // delete X coef
      //should find row in uniqueS that corresponds to person i
      ind_dummy = find(uniqueS.col(0)==Sy(i) && uniqueS.col(1)==Sx(i));

      xPiPars.shed_row(ind_dummy(0)); // CHECK --- NEEDED TO DO THIS BECAUSE ind_dummy is uvec
      xMuPars.shed_row(ind_dummy(0));
      xSigPars.shed_row(ind_dummy(0));

      //relabel X cluster
      for(int j = 0; j < Sx.size(); j++) {
	if(Sy(j) == Sy(i) && Sx(j) > Sx(i))
	  Sx(j) = Sx(j) - 1;
      }

      for(int j = 0; j < uniqueS.n_rows; j++) {
	if(uniqueS(j,0) == Sy(i) && uniqueS(j,1) > Sx(i))
	  uniqueS(j,1) = uniqueS(j,1) - 1;
      }

      //relabel Y cluster (if needed)
      if(num_dummy2==1) {
	for(int j = 0; j < Sy.size(); j++) {
	  if(Sy(j) > Sy(i))
	    Sy(j) = Sy(j) - 1;
	}

	for(int j = 0; j < uniqueS.n_rows; j++) {
	  if(uniqueS(j,0) > Sy(i))
	    uniqueS(j,0) = uniqueS(j,0) - 1;
	}
      }

      uniqueS.shed_row(ind_dummy(0)); //get rid of row
    }

    // NEED TO DELETE ROW OF Sy and Sx??? Think so
    Sy.shed_row(i); Sx.shed_row(i);


    // recalculate number of unique clusters
    numY = betaY.n_rows;
    numTotalCluster = xMuPars.n_rows;
    int totalposs = numY+numTotalCluster+1;
    //Rprintf("Y clusters, total clusters, Total possible clusters: %d,%d,%d\n",numY,numTotalCluster,totalposs);
    vec probs(totalposs);
    count=0;

    int njwoi, nljwoi; //counts for # in appropriate Y and X cluster,
                       //excluding the ith person

    int numXj; //# of X clusters within each Y cluster


    for(int j = 0; j < numY; j++) {
      //FILL IN PROBS FOR EXISTING CLUSTERS

      //get count of number of X clusters within jth Y cluster
      ind_dummy = find(uniqueS.col(0) == (j+1));
      numXj = ind_dummy.size();
      //Rprintf("NumX cluster within Y cluster: %d\n",numXj);

      //get number of subjects within jth cluster
      ind_dummy = find(Sy == (j+1));
      njwoi = ind_dummy.size();

      //Rprintf("njwoi: %d\n",njwoi);
      // likelihood for each existing Y cluster
      //likeregy = R::dbinom(y(i), 1 , expit( dot( matX.row(i),betaY.row(j) ) ) , 0);
      likeregy = R::dnorm(y(i),dot(matX.row(i),betaY.row(j)),sqrt(sig2(j)),0);

      //Rprintf("likelihood for Y: %.2f\n",likeregy);
      for(int k = 0; k < numXj; k++) {
	prodx = 1;
	prodx2 = 1;

	ind_dummy = find(Sy == (j+1) && Sx == (k+1));
	nljwoi = ind_dummy.size();
	//Rprintf("nljwoi: %d\n",nljwoi);

	//likelihood for binary covariates
	for(int l = 0; l < ptx+p1; l++) {
	  prodx = prodx*R::dbinom(X(i,l), 1, xPiPars(count, l), 0);
	}

	//Rprintf("Prodx: %.2f\n",prodx);

	//likelihood for continuous covariates
	for(int l = 0; l < p2; l++) {
	  prodx2 = prodx2*R::dnorm(X(i,ptx+p1+l), xMuPars(count,l), sqrt(xSigPars(count,l)), 0 );
	}

	//Rprintf("Prodx2: %.2f\n",prodx2);

	probs(count) = ((njwoi*nljwoi)/(njwoi+alphapsi))*likeregy*prodx*prodx2;
	//Rprintf("\t current count and prob: %d, %.2f\n",count,probs(count));
	count++;
      }
    }


    for(int j = 0; j < numY; j++) {
      //FILL IN PROBS FOR NEW X CLUSTERS IN EXISTING Y CLUSTERS

      ind_dummy = find(Sy == (j+1));
      njwoi = ind_dummy.size();

      // likelihood for each existing Y cluster
      //likeregy = R::dbinom(y(j), 1 , expit( dot( matX.row(i),betaY.row(j) ) ) , 0);
      // NEW FOR CONTINUOUS
      likeregy = R::dnorm(y(i),dot(matX.row(i),betaY.row(j)),sqrt(sig2(j)),0);
      probs(numTotalCluster + j) = ( (njwoi*alphapsi) / (njwoi+alphapsi) ) * likeregy*h0i(i);
      //Rprintf("New X cluster\n");
      //Rprintf("\t current count and prob: %d, %.2f\n",numTotalCluster+j,probs(numTotalCluster+j));
    }

    probs(numY+numTotalCluster) = alphatheta*h0y(i)*h0i(i); //prob for new Y cluster
    //Rprintf("New Y cluster prob\n");
    //Rprintf("\t current count and prob: %d, %.2f\n",numY+numTotalCluster,probs(numY+numTotalCluster));

    //USE MULTINOMIAL DISTRIBUTION TO CHOOSE NEW CLUSTER
    newCluster = rmultinomF(probs);
    probs.zeros();

    //Rprintf("The new cluster is: %d\n", newCluster);
    //need to map this integer to one of the clusters (or a new cluster)

    if(newCluster<=numTotalCluster)
      {
	Sy.insert_rows(i,1);
	Sy(i) = uniqueS(newCluster-1,0);
	Sx.insert_rows(i,1);
	Sx(i) = uniqueS(newCluster-1,1);
      }
    else
      {
	//find out whether this is a new Y cluster or X cluster
	if(newCluster == (numTotalCluster+numY+1))
	  {

	    Sx.insert_rows(i,1);
	    Sx(i) = 1;

	    // need functions rmvn, newbetafunction, updatevar, updatemean

	    //betadraw = mvrnorm(betainit, diagbetacov0);
	    //newbeta = newbetafunction(betadraw,trans(matX.row(i)),y(i),betainit,diagbetacov0);
	    newsig = newregvar(trans(matX.row(i)),y(i),betaa0,betab0,betainit);
	    newbeta = newbet(trans(matX.row(i)),y(i),newsig,betainit);
	    betaY.insert_rows(numY,newbeta.t());
	    sig2.insert_rows(numY,1);
	    sig2(numY) = newsig;

	    for(int j = 0; j<(p1+ptx); j++) {
	      newpipars(j) = R::rbeta(X(i,j)+a0,1-X(i,j)+b0);
	    }

	    for(int j = 0; j < p2; j++) {
	      newsigpars(j) = updatevar(X(i,(p1+ptx+j)),nu0,tau,c0,mu0);
	      newmupars(j) = updatemean(X(i,(p1+ptx+j)),tau,c0,mu0);
	    }
	    xPiPars.insert_rows(numTotalCluster,newpipars);
	    xSigPars.insert_rows(numTotalCluster,newsigpars);
	    xMuPars.insert_rows(numTotalCluster,newmupars);


	    Sy.insert_rows(i,1);
	    Sy(i) = Sy.max()+1;
	    uniqueS.insert_rows(numTotalCluster,1);
	    uniqueS(numTotalCluster,0) = Sy(i);
	    uniqueS(numTotalCluster,1) = Sx(i);
	  }
	else
	  {
	    //Y cluster should be numCluster - numTotalCluster
	    //X cluster should be max +1 of X clusters within Y

	    Sy.insert_rows(i,1);
	    Sy(i) = newCluster - numTotalCluster;
	    ind_dummy = find(uniqueS.col(0) == Sy(i));
	    Sx.insert_rows(i,1);
	    Sx(i) = ind_dummy.size() + 1;

	    // update parameters
	    for(int j = 0; j<(p1+ptx); j++) {
	      newpipars(j) = R::rbeta(X(i,j)+a0,1-X(i,j)+b0);
	    }

	    for(int j = 0; j < p2; j++) {
	      newsigpars(j) = updatevar(X(i,(p1+ptx+j)),nu0,tau,c0,mu0);
	      newmupars(j) = updatemean(X(i,(p1+ptx+j)),tau,c0,mu0);
	    }

	    ind_dummy = find(uniqueS.col(0) <= Sy(i));
	    num_dummy = ind_dummy.size();

	    xPiPars.insert_rows(num_dummy,newpipars);
	    xSigPars.insert_rows(num_dummy,newsigpars);
	    xMuPars.insert_rows(num_dummy,newmupars);

	    uniqueS.insert_rows(num_dummy,1);
	    uniqueS(num_dummy,0) = Sy(i);
	    uniqueS(num_dummy,1) = Sx(i);
	  }
      }


    // as necessary add back into Sy, Sx, uniqueS, etc...
    //Rprintf("New cluster for person i: %d, %d\n\n",Sy(i),Sx(i));
    //Rcout << uniqueS << std::endl;
    //Rcout << betaY << std::endl;
    //Rprintf("Max Sy, Max uniqueS: %d,%d\n",Sy.max(),uniqueS.col(0).max());
  }



  return List::create(_["Sy"]  = Sy,
		      _["Sx"] = Sx,
		      _["uniqueS"] = uniqueS,
		      _["beta"]    = betaY,
		      _["sig2"]    = sig2,
		      _["pipars"]  = xPiPars,
		      _["mupars"]  = xMuPars,
		      _["sigpars"] = xSigPars);
}
