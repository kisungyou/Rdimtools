#ifndef _Rdimtools_METHODS_LINEAR_H
#define _Rdimtools_METHODS_LINEAR_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <stdlib.h>

using namespace Rcpp;
using namespace arma;

/*
 * MAIN METHODS
 * 01. PCA
 * 02. MDS
 * 03. MDS given D
 * 04. ICA
 * 05. RNDPROJ
 * 06. FA
 * 07.
 * 08. NPE
 * 09. OLPP
 * 10. BPCA
 * 11. EXTLPP
 * 12. LSPP
 * 13. KMMC
 * 14. LFDA
 * 15. NNPROJMAX & NNPROJMIN
 * 16. NNEMBEDMIN
 */

Rcpp::List method_pca(arma::mat& psdX);                             // 01. PCA
Rcpp::List method_mds(arma::mat& centerX);                          // 02. MDS
Rcpp::List method_mdsD(arma::mat& D);                               // 03. MDS given D
Rcpp::List method_ica(arma::mat& X, const int C,
                      const int maxiter, const double tol,
                      const int tnum, const double tpar, bool sym); // 04. ICA
Rcpp::List method_rpgauss(arma::mat& X, const int k);               // 05. RNDPROJ
// Rcpp::List method_fa(arma::mat& X, const int k, const int maxiter,  // 06. FA
                     // const double tolerance);
Rcpp::List method_npe(arma::mat& X, arma::mat& W);                  // 08. NPE
arma::mat method_olpp(arma::mat& X, arma::mat& S, const int ndim);  // 09. OLPP
arma::mat auxiliary_outer(arma::colvec x, arma::colvec y);          // 10. BPCA
Rcpp::List method_bpca(arma::mat& T, const double reltol, const int maxiter);
arma::mat method_trfextlpp(arma::mat& D, double a, double b);       // 11. EXTLPP
arma::mat method_lspp_computeW(arma::mat& S, arma::vec& svec);      // 12. LSPP
arma::vec method_kmmcvec(arma::mat& X, arma::mat& partmat,          // 13. KMMC
                         double param);
double method_lfda_maximaldistance(arma::rowvec& tvec,              // 14. LFDA
                                   arma::mat& tmat);
arma::mat method_nnprojmax(arma::mat& C, arma::mat& Uinit,          // 15. NNPROJMAX & NNPROJMIN
                      const double tol, const int maxiter);
arma::mat method_nnprojmin(arma::mat& C, arma::mat& Uinit, const double tol, const int maxiter);
arma::mat method_nnembedmin(arma::mat& M, arma::mat& Yinit,         // 16. NNEMBEDMIN
                            const double tol, const int maxiter);
arma::vec method_spufs(arma::mat& X, arma::mat Ls, double alpha, double beta, double epsilon);
arma::vec method_lspe(arma::mat X, const int d, double alpha, double beta, arma::mat L);
arma::vec method_disr(arma::mat& D, double lbd1, double lbd2);
arma::vec method_rsr(arma::mat X, double lbd, double verysmall);
arma::vec method_nrsr(arma::mat X, double lbd, double verysmall, double p);
arma::vec method_scoresum(arma::mat &X, arma::mat &S);
#endif
