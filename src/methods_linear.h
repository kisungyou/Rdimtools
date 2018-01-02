#ifndef _Rdimtools_METHODS_LINEAR_H
#define _Rdimtools_METHODS_LINEAR_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <stdlib.h>

using namespace Rcpp;
using namespace arma;

// TO DO : method_mdsSD - sparse distance matrix !
// 01. PCA
Rcpp::List method_pca(arma::mat& psdX);
// 02. MDS
Rcpp::List method_mds(arma::mat& centerX);
// 03. MDS given D
Rcpp::List method_mdsD(arma::mat& D);
// 04. ICA
Rcpp::List method_ica(arma::mat& X, const int C, const int maxiter, const double tol, const int tnum, const double tpar, bool sym);
// 05. RNDPROJ
Rcpp::List method_rpgauss(arma::mat& X, const int k);
// 06. FA
Rcpp::List method_fa(arma::mat& X, const int k, const int maxiter, const double tolerance);
// 07. LPP
Rcpp::List method_lpp(arma::mat& X, arma::mat& W);
// 08. NPE
Rcpp::List method_npe(arma::mat& X, arma::mat& W);
// 09. OLPP
arma::mat method_olpp(arma::mat& X, arma::mat& S, const int ndim);
// 10. BPCA
arma::mat auxiliary_outer(arma::colvec x, arma::colvec y);
Rcpp::List method_bpca(arma::mat& T, const double reltol, const int maxiter);
// 11. EXTLPP
arma::mat method_trfextlpp(arma::mat& D, double a, double b);
// 12. LSPP
arma::mat method_lspp_computeW(arma::mat& S, arma::vec& svec);
// 13. KMMC
arma::vec method_kmmcvec(arma::mat& X, arma::mat& partmat, double param);
// 14. LFDA
double method_lfda_maximaldistance(arma::rowvec& tvec, arma::mat& tmat);

#endif
