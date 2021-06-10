#include <RcppArmadillo.h>
#include "ver1_computation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// 01. PCA
// 02. MDS
// 03. FA
// 04. LMDS
// 05. SPCA

// 01. PCA ====================================================================
// [[Rcpp::export]]
Rcpp::List dt_pca(arma::mat& X, int ndim, bool cor){
  // preliminary --------------------------------------------------------------
  // parameters
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.pca : 'ndim' should be in [1,ncol(X)).");
  }

  // main computation ---------------------------------------------------------
  // 1. construct
  arma::mat psdX(P,P);
  if (cor==true){
    psdX = arma::cor(X);
  } else {
    psdX = arma::cov(X);
  }
  // 2. eigendecomposition
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, psdX); // ascending order issue

  arma::vec vars = arma::reverse(eigval.tail(ndim));
  arma::mat proj = arma::fliplr(eigvec.tail_cols(ndim));

  // 3. computation
  arma::mat Y = X*proj;

  // repot ---------------------------------------------------------------------
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("vars") = vars,
      Rcpp::Named("projection") = proj,
      Rcpp::Named("algorithm") = "linear:PCA"
  ));
}

// 02. MDS ================================================================
// [[Rcpp::export]]
Rcpp::List dt_mds(arma::mat& X, int ndim){
  // preliminary --------------------------------------------------------------
  // parameters
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.mds : 'ndim' should be in [1,ncol(X)).");
  }
  arma::mat Xt = arma::trans(X);

  // main computation ---------------------------------------------------------
  // 1. eigendecomposition
  arma::mat U, V; arma::vec s;
  arma::svd(U,s,V,Xt);
  // 2. compute Y first
  arma::mat Y = V.head_cols(ndim)*arma::diagmat(s.head(ndim));
  // 3. compute Projection
  arma::mat LHS  = Xt*X;
  arma::mat RHS  = Xt*Y;
  arma::mat proj = v2aux_solproj(LHS,RHS);

  // wrap and report ---------------------------------------------------
  // trfinfo
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("projection") = proj,
      Rcpp::Named("algorithm") = "linear:MDS"
  ));
}

// 03. FA ======================================================================
// [[Rcpp::export]]
Rcpp::List dt_fa(arma::mat& X, int ndim, int maxiter, double tolerance){
  // preliminary --------------------------------------------------------------
  // parameters
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.fa : 'ndim' should be in [1,ncol(X)).");
  }
  arma::mat Xt = arma::trans(X);

  // main computation ---------------------------------------------------------
  // 1. pass onto v2aux_fa
  Rcpp::List output = v2aux_fa(Xt,ndim,maxiter,tolerance);

  // 2. after works
  arma::mat Y    = output["Z"];
  arma::mat LHS  = Xt*X;
  arma::mat RHS  = Xt*Y.t();
  arma::mat projtmp = arma::solve(LHS,RHS);
  arma::mat proj    = v2aux_adjproj(projtmp);

  arma::mat outL = output["L"];
  arma::vec outP = output["Pvec"];


  // wrap and report ---------------------------------------------------
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y.t(),
      Rcpp::Named("noise")   = outP,
      Rcpp::Named("projection") = proj,
      Rcpp::Named("loadings")   = outL,
      Rcpp::Named("algorithm") = "linear:FA"
  ));
}

// 04. LMDS -===================================================================
// [[Rcpp::export]]
Rcpp::List dt_lmds(arma::mat& X, int ndim, int npts){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.lmds : 'ndim' should be in [1,ncol(X)).");
  }
  if ((npts<2)||(npts>N)){
    throw std::invalid_argument("* do.lmds : the number of landmark points is not valid.");
  }

  // main computation ---------------------------------------------------------
  // 1. random selection
  arma::uvec idselect = arma::randperm(N,npts);
  // 2. do the pca
  arma::mat subX = X.rows(idselect);
  arma::mat psdX = arma::cov(subX);
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, psdX); // ascending order issue
  arma::vec vars = arma::reverse(eigval.tail(ndim));
  arma::mat proj = arma::fliplr(eigvec.tail_cols(ndim));

  // 3. computation
  arma::mat Y = X*proj;

  // wrap and report ---------------------------------------------------
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("projection") = proj,
      Rcpp::Named("algorithm") = "linear:LMDS"
  ));
}



// 05. SPCA ====================================================================
// auxiliary functions for dt_spca
arma::vec dt_spca_rk1vec(arma::mat& X){
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, X);

  arma::vec output = eigvec.tail_cols(1);
  return(output);
}
arma::mat dt_spca_deflation(arma::mat& Sig, arma::vec& Vec){
  arma::mat term1 = Sig*(Vec*Vec.t())*Sig; // column-vector convention
  double term2 = arma::dot((Sig*Vec), Vec);
  arma::mat output = Sig - (term1/term2);
  return(output);
}
// main function for spca
// [[Rcpp::export]]
Rcpp::List dt_spca(arma::mat& X, int ndim, double mu, double rho, const double abstol, const double reltol, const int maxiter){
  // preliminary --------------------------------------------------------------
  // parameters
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.spca : 'ndim' should be in [1,ncol(X)).");
  }
  if (mu < arma::datum::eps){
    throw std::invalid_argument("* do.spca : 'mu' should be a positive real number.");
  }
  if (rho < arma::datum::eps){
    throw std::invalid_argument("* do.spca : 'rho' should be a positive real number.");
  }

  // main computation ---------------------------------------------------------
  // 1. prepare for name changes in ADMM's SPCA function.
  int numpc = ndim;
  arma::mat Sigma = arma::cov(X);
  arma::mat basis(P,numpc,fill::zeros);
  arma::mat tmpX;
  arma::vec solvec;

  // 2. iteration
  for (int i=0;i<numpc;i++){
    // 2-1. compute tmpX
    tmpX   = v2aux_spca(Sigma, reltol, abstol, maxiter, mu, rho);
    // 2-2. compute solvec
    solvec = dt_spca_rk1vec(tmpX); // possibly type conversion required
    basis.col(i) = solvec;
    // 2-3. update Sigma
    Sigma = dt_spca_deflation(Sigma, solvec);
  }

  // 3.   compute basis & Y
  arma::mat proj = v2aux_adjproj(basis);
  arma::mat Y    = X*proj;

  // wrap and report ---------------------------------------------------
  // trfinfo
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("projection") = proj,
      Rcpp::Named("algorithm") = "linear:SPCA"
  ));
}




