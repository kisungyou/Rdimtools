#include <RcppArmadillo.h>
#include "ver2_classes.h"
#include "ver2_computation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/*
 * 01. dt_pca : PCA
 * 02. dt_fa  : FA
 */



// 01. PCA ====================================================================
// [[Rcpp::export]]
Rcpp::List dt_pca(const arma::mat& X, int ndim, std::string ptype, bool cor){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.pca : 'ndim' should be in [1,ncol(X)).");
  }

  // preprocessing
  ClPreproc init(ptype);
  arma::mat premat    = init.MainFunc(X);
  arma::rowvec mymean = premat.row(0);
  arma::mat    mymult = premat.rows(1,P);
  std::string  mytype = init.GetType();

  // data prep
  arma::mat Xtmp(N,P,fill::zeros);
  for (int n=0;n<N;n++){
    Xtmp.row(n) = X.row(n) - mymean;
  }
  Xtmp = Xtmp*mymult;

  // main computation ---------------------------------------------------------
  // 1. construct
  arma::mat psdX(P,P);
  if (cor==true){
    psdX = arma::cor(Xtmp);
  } else {
    psdX = arma::cov(Xtmp);
  }
  // 2. eigendecomposition
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, psdX); // ascending order issue

  arma::vec vars = arma::reverse(eigval.tail(ndim));
  arma::mat proj = arma::fliplr(eigvec.tail_cols(ndim));

  // 3. computation
  arma::mat Y = Xtmp*proj;

  // wrap and report ---------------------------------------------------
  // trfinfo
  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("type")=mytype,
    Rcpp::Named("mean")=mymean,
    Rcpp::Named("multiplier")=mymult,
    Rcpp::Named("algtype")="linear"
  );
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("vars") = vars,
      Rcpp::Named("projection") = proj,
      Rcpp::Named("trfinfo") = info
  ));
}

// 02. FA =====================================================================
// [[Rcpp::export]]
Rcpp::List dt_fa(const arma::mat& X, int ndim, std::string ptype, int maxiter, double tolerance){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.fa : 'ndim' should be in [1,ncol(X)).");
  }

  // preprocessing
  ClPreproc init(ptype);
  arma::mat premat    = init.MainFunc(X);
  arma::rowvec mymean = premat.row(0);
  arma::mat    mymult = premat.rows(1,P);
  std::string  mytype = init.GetType();

  // data prep
  arma::mat pX(N,P,fill::zeros);
  for (int n=0;n<N;n++){
    pX.row(n) = X.row(n) - mymean;
  }
  pX = pX*mymult;

  // main computation ---------------------------------------------------------
  // 1. pass onto v2aux_fa
  Rcpp::List output = v2aux_fa(arma::trans(pX),ndim,maxiter,tolerance);

  // 2. after works
  arma::mat Y = output["Z"];

  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("type")=mytype,
    Rcpp::Named("mean")=mymean,
    Rcpp::Named("multiplier")=mymult,
    Rcpp::Named("algtype")="linear"
  );

  arma::mat LHS  = arma::trans(pX)*pX;
  arma::mat RHS  = arma::trans(pX)*Y.t();
  arma::mat projtmp = arma::solve(LHS,RHS);
  arma::mat proj    = v2aux_adjproj(projtmp);

  arma::mat outL = output["L"];
  arma::vec outP = output["Pvec"];


  // wrap and report ---------------------------------------------------
  // trfinfo
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y.t(),
      Rcpp::Named("trfinfo") = info,
      Rcpp::Named("noise")   = outP,
      Rcpp::Named("projection") = proj,
      Rcpp::Named("loadings")   = outL
  ));
}
