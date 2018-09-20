#ifndef _Rdimtools_METHODS_NONLINEAR_H
#define _Rdimtools_METHODS_NONLINEAR_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// 1. SNE
arma::mat method_sne(arma::mat& P, const int ndim, const double eta,
                     const int maxiter, double jitter, double decay,
                     const double momentum);
// 2. Symmetric SNE
arma::mat method_snesym(arma::mat& P, const int ndim, const double eta,
                     const int maxiter, double jitter, double decay,
                     const double momentum);
// 3. t-SNE
arma::mat method_tsne(arma::mat& P, const int ndim, const double eta,
                        const int maxiter, double jitter, double decay,
                        const double momentum);
// 4. Eigenmaps
Rcpp::List method_eigenmaps(arma::mat& W);
// 5. Sammon Mapping
arma::mat method_sammon(arma::mat& X, arma::mat& init);
// 6. LLE
arma::vec method_lleW(arma::mat& mat_tgt, arma::vec& vec_tgt, const double regparam);
// 7. LLE with automatic choice
Rcpp::List method_lleWauto(arma::mat& mat_tgt, arma::vec& vec_tgt);
// 8. LLE M
Rcpp::List method_lleM(arma::mat& W);
// 9. REE
double method_ree_cost(arma::mat W, arma::mat D, arma::mat B);
arma::mat method_ree_subgradient(arma::mat B, arma::mat W, arma::mat D);
Rcpp::List method_ree(arma::mat& B, arma::mat& W, arma::mat& D, const double initc,
                      const double abstol, const int maxiter);
// 10. SPE
arma::mat method_spe(arma::mat& R, arma::mat& iX, const int C, const int S, double lambda, double drate, arma::mat matselector);
arma::mat method_ispe(arma::mat& R, arma::mat& iX, const int C, const int S,
                      double lambda, double drate, arma::mat matselector, const double cutoff);
// 11. CRCA : Curvilinear Component Analysis
Rcpp::List method_crca(arma::mat& Xij, arma::mat& Yinit, double lambda, double alpha, const int maxiter, const double tolerance, arma::vec& vecselector);

#endif
