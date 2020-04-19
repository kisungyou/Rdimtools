#ifndef _Rdimtools_VER2_COMPUTATION_H
#define _Rdimtools_VER2_COMPUTATION_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

/* Part 1. Methods Related
 *    v2aux_fa     : Factor Analysis
 *    v2aux_pca    : Principal Component Analysis
 *
 * Part 2. Generic
 *    v2aux_pdist   : pairwise distance in L2 norm
 *    v2aux_adjproj : replacing 'aux.adjprojection'
 *
 */

// Part 1. Methods Related
Rcpp::List v2aux_fa(arma::mat& X, const int k, const int maxiter, const double tolerance);
arma::mat  v2aux_pca(arma::mat& X, int ndim);

// Part 2. Generic
arma::mat  v2aux_pdist(arma::mat& X);
arma::mat  v2aux_adjproj(arma::mat& X);


#endif
