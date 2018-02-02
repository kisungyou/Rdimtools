#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List aux_mimick_geigen(arma::mat& A, arma::mat& B){
  arma::cx_vec eigval;
  arma::cx_mat eigmat;

  eig_pair(eigval, eigmat, A, B);
  return Rcpp::List::create(Rcpp::Named("values")=eigval,
                            Rcpp::Named("vectors")=eigmat);
}
