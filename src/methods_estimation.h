#ifndef _Rdimtools_METHODS_ESTIMATION_H
#define _Rdimtools_METHODS_ESTIMATION_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// 1. boxcounting
arma::mat methods_boxcount(arma::mat& X, arma::vec& Imin, const double currentr);
// 2. numerical derivative
arma::vec aux_numderiv(arma::vec& x, arma::vec& y);

#endif
