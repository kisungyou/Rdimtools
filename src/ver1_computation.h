#ifndef _Rdimtools_VER1_COMPUTATION_H
#define _Rdimtools_VER1_COMPUTATION_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

/* Part 1. Methods Related
 *    v2aux_fa     : FA
 *    v2aux_pca    : PCA
 *    v2aux_spca   : Sparse PCA
 *
 * Part 2. Generic
 *    v2aux_pdist    : pairwise distance in L2 norm
 *    v2aux_adjproj  : replacing 'aux.adjprojection'
 *    v2aux_solproj  : solve linear system and return projection with QR
 *    v2aux_knn      : return the index of k-nearest neighbors
 *    v2aux_fid2proj : replicate 'aux.featureindicator' function
 *    v2aux_pagerank : page rank algorithm
 */

// Part 1. Methods Related
Rcpp::List v2aux_fa(arma::mat& X, const int k, const int maxiter, const double tolerance);
arma::mat  v2aux_pca(arma::mat& X, int ndim);
arma::mat  v2aux_spca(arma::mat& Sigma, const double reltol, const double abstol, const int maxiter, double mu, double rho);

// Part 2. Generic
arma::mat  v2aux_pdist(arma::mat& X);
arma::mat  v2aux_adjproj(arma::mat& X);
arma::mat  v2aux_solproj(arma::mat& LHS, arma::mat& RHS);
arma::umat v2aux_knn(arma::mat& X, int k);
arma::mat  v2aux_fid2proj(int p, int ndim, arma::uvec idxvec);
arma::vec  v2aux_pagerank(arma::mat& A);

#endif
