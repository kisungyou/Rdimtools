#ifndef _Rdimtools_METHODS_NONLINEAR_H
#define _Rdimtools_METHODS_NONLINEAR_H

#define ARMA_NO_DEBUG

#include <RcppDist.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. SNE
arma::mat method_sne(arma::mat& P, int ndim0, double eta0,
                     int maxiter0, double jitter0, double decay0,
                     double momentum0);
// 2. Symmetric SNE
arma::mat method_snesym(arma::mat& P, int ndim0, double eta0,
                        int maxiter0, double jitter0, double decay0,
                        double momentum0);
// 3. t-SNE
arma::mat method_tsne(arma::mat& P, int ndim0, double eta0,
                      int maxiter0, double jitter0, double decay0,
                      double momentum0);
// 4. Eigenmaps
Rcpp::List method_eigenmaps(arma::mat& W);
// 5. Sammon Mapping - use MASS
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

// 12. BMDS
double bmds_compute_SSR(arma::mat &D, arma::mat &Delta);
double bmds_compute_SSR_xmat(arma::mat &D, arma::mat &Xnew);
arma::mat bmds_compute_pdmat(arma::mat &X);
arma::mat bmds_crotX(arma::mat X);
arma::rowvec bmds_update_xvec(arma::mat D, arma::mat X, int id, double sigma2, double constant, arma::mat Lbdmat);
double my_invgamma(double alpha, double beta);
double my_dinvgamma(double x, double alpha, double beta);
Rcpp::List main_bmds(arma::mat D, arma::mat X0, double sigg0,
                     double a, double alpha, int maxiter, double constant, bool verbose,
                     arma::vec betas);
#endif
