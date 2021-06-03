#include <RcppArmadillo.h>
#include "ver1_computation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// 01. MMDS

// 01. MMDS : metric MDS -------------------------------------------------------
arma::mat dt_mmds_core(arma::mat D, int ndim, int maxiter, double abstol){
  // initialization with CMDS
  int N = D.n_rows;
  arma::mat D2 = arma::pow(D, 2.0);
  arma::mat J  = arma::eye<arma::mat>(N,N) - (arma::ones<arma::mat>(N,N)/(static_cast<double>(N)));
  arma::mat B  = -0.5*J*D2*J;  arma::vec eigval;  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, B);

  arma::mat old_X = eigvec.tail_cols(ndim)*arma::diagmat(arma::sqrt(eigval.tail(ndim)));
  arma::mat old_D(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      old_D(i,j) = arma::norm(old_X.row(i)-old_X.row(j),2);
      old_D(j,i) = old_D(i,j);
    }
  }

  arma::mat new_X(N,ndim,fill::zeros);
  arma::mat new_D(N,N,fill::zeros);
  double old_cost = 0.0;
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      old_cost += std::pow(D(i,j)-old_D(i,j), 2.0);
    }
  }
  double new_cost = 0.0;

  arma::mat BZ(N,N,fill::zeros);
  double dijZ   = 0.0;
  double epsthr = 100*arma::datum::eps;
  double inctol = 0.0;
  double diagsum = 0.0;

  for (int it=0; it<maxiter; it++){
    // compute updater B(Z); (Borg p155); Z=old_X
    BZ.fill(0.0);
    // off-diagonal first
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        // dijZ = arma::norm(old_X.row(i)-old_X.row(j),2);
        dijZ = old_D(i,j);
        if (dijZ > epsthr){
          BZ(i,j) = -D(i,j)/dijZ;
          BZ(j,i) = BZ(i,j);
        }
      }
    }
    // diagonal part
    for (int i=0; i<N; i++){
      diagsum = -arma::accu(BZ.row(i));
      BZ(i,i) = diagsum;
    }
    // updater
    new_X    = (BZ*old_X)/(static_cast<double>(N));
    new_D.fill(0.0);
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        new_D(i,j) = arma::norm(new_X.row(i)-new_X.row(j),2);
        new_D(j,i) = new_D(i,j);
      }
    }
    new_cost = 0.0;
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        new_cost += std::pow(D(i,j)-new_D(i,j), 2.0);
      }
    }

    inctol   = old_cost - new_cost;
    old_X    = new_X;
    old_cost = new_cost;
    if (inctol < abstol){
      break;
    }
  }
  return(old_X);
}
// [[Rcpp::export]]
Rcpp::List dt_mmds(arma::mat& X, int ndim, int maxiter, double abstol){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= X.n_cols)){
    throw std::invalid_argument("* do.mmds : 'ndim' should be in [1,ncol(X)).");
  }

  // computation ---------------------------------------------------------------
  // distance computation
  arma::mat D(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      D(i,j) = arma::norm(X.row(i)-X.row(j),2);
      D(j,i) = D(i,j);
    }
  }
  // compute
  arma::mat embed = dt_mmds_core(D, ndim, maxiter, abstol);

  // wrap and report ---------------------------------------------------
  return(Rcpp::List::create(
      Rcpp::Named("Y") = embed,
      Rcpp::Named("algorithm") = "nonlinear:mmds"
  ));
}
