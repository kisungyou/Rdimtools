#ifndef _Rdimtools_METHODS_HADNYTOOLS_H
#define _Rdimtools_METHODS_HADNYTOOLS_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <stdlib.h>

using namespace Rcpp;
using namespace arma;


arma::mat handy_plus(arma::mat& X);                                    // 1. handy_plus         : from X, find X+
arma::mat handy_hadamartABC(arma::mat& A, arma::mat& B, arma::mat& C); // 2. handy_hadamardABC  : A*B/C
arma::mat handy_hadamartABCsqrt(arma::mat& A, arma::mat& B,            //    handy_hadamardABCsqrt
                                arma::mat& C);

#endif
