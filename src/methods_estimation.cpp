#include <RcppArmadillo.h>
#include "methods_estimation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


/*
 * 1. boxcounting
 *    generate integer-valued data labeling
 */
// [[Rcpp::export]]
arma::mat methods_boxcount(arma::mat& tX, arma::vec& Imin, double currentr){
  // 1-1. basic settings
  int d  = tX.n_rows;
  int n  = tX.n_cols;
  int dI = Imin.n_elem;
  if (d != dI){
    Rcpp::stop("ERROR : dimension not matching.");
  }

  // 1-2. main iteration
  vec target;
  mat intmat(n,d,fill::zeros);
  for (int i=0;i<n;i++){
    // 1-2-1. select target data
    target = tX.col(i) - Imin;
    // 1-2-2. find integer labeling
    for (int j=0;j<d;j++){
      intmat(i,j) = floor(as_scalar(target(j)/currentr));
    }
  }

  // 1-3. return output
  return(intmat);
}

/*
 * 2. numderiv : use all forward difference except the last one
 */
// [[Rcpp::export]]
arma::vec aux_numderiv(arma::vec& x, arma::vec& y){
  int n = x.n_elem;
  arma::vec deriv(n,fill::zeros);

  deriv(n-1) = (y(n-1)-y(n-2))/(x(n-1)-x(n-2));
  for (int i=0;i<(n-1);i++){
    deriv(i) = (y(i+1)-y(i))/(x(i+1)-x(i));
  }
  return(deriv);
}
