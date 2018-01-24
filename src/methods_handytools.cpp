#include <RcppArmadillo.h>
#include "methods_handytools.h"
#include <stdlib.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


/*
 * 1. handy_plus : from X, find X+
 */
// [[Rcpp::export]]
arma::mat handy_plus(arma::mat& X){
  const int n = X.n_rows;
  const int p = X.n_cols;

  arma::mat output(n,p,fill::zeros);
  for (int i=0;i<n;i++){
    for (int j=0;j<p;j++){
      if (X(i,j)>=0){
        output(i,j) = X(i,j);
      }
    }
  }
  return(output);
}

/*
 * 2. handy_hadamardABC      : A*B/C       in R notation
 *    handy_hadamardABCsqrt  : A*sqrt(B/C) in R notation
 */
// [[Rcpp::export]]
arma::mat handy_hadamartABC(arma::mat& A, arma::mat& B, arma::mat& C){
  const int n = A.n_rows;
  const int p = A.n_cols;

  arma::mat output(n,p,fill::zeros);
  for (int i=0;i<n;i++){
    for (int j=0;j<p;j++){
      if (C(i,j)!=0){
        output(i,j) = A(i,j)*B(i,j)/C(i,j);
      }
    }
  }
  return(output);
}
// [[Rcpp::export]]
arma::mat handy_hadamartABCsqrt(arma::mat& A, arma::mat& B, arma::mat& C){
  const int n = A.n_rows;
  const int p = A.n_cols;

  arma::mat output(n,p,fill::zeros);
  for (int i=0;i<n;i++){
    for (int j=0;j<p;j++){
      if (C(i,j)!=0){
        output(i,j) = A(i,j)*(static_cast<double>(std::sqrt(B(i,j)/C(i,j))));
      }
    }
  }
  return(output);
}
