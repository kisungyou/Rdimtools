#include <RcppArmadillo.h>
#include "ver2_computation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat v2aux_pca(const arma::mat& X, int ndim){
  int N = X.n_rows;
  int P = X.n_cols;

  arma::mat psdX = arma::cov(X);
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, psdX); // ascending order issue

  arma::rowvec xmean = arma::mean(X, 0);
  arma::mat proj = arma::fliplr(eigvec.tail_cols(ndim));

  arma::mat output(N,ndim,arma::fill::zeros);
  for (int n=0;n<N;n++){
    output.row(n) = (X.row(n)-xmean)*proj;
  }
  return(output);
}
arma::mat v2aux_pdist(const arma::mat& X){
  int N = X.n_rows;
  int P = X.n_cols;

  arma::mat output(N,N,fill::zeros);
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      output(i,j) = arma::norm(X.row(i)-X.row(j),2);
      output(j,i) = output(i,j);
    }
  }
  return(output);
}
Rcpp::List v2aux_fa(const arma::mat& X, const int k, const int maxiter, const double tolerance){
  // 6-1. basic settings
  const int p = X.n_rows;
  const int n = X.n_cols;

  // 6-2. other defined materials
  mat Ez(k,n);
  mat beta(k,p);
  mat LPhi(p,p);
  mat eyeK(k,k,fill::eye);
  mat Pinv(p,p,fill::zeros);
  mat EzzSum(k,k);
  double inctol = 0;

  // 6-3. initialization
  arma::mat Lold(p,k,fill::randn);
  arma::mat Lnew(p,k,fill::zeros);
  arma::vec Pold(p,fill::ones);
  arma::mat Ptmp(p,p,fill::zeros);
  arma::vec Pnew(p,fill::zeros);

  // 6-4. main iteration
  for (int it=0;it<maxiter;it++){
    // 6-4-E1. common again
    Pinv = diagmat(1/Pold);
    // 6-4-E2. inverse of LL^T+\Phi : LPhi
    LPhi = Pinv - Pinv*Lold*solve(eyeK+Lold.t()*Pinv*Lold,Lold.t()*Pinv);
    // 6-4-E3. beta
    beta = Lold.t()*LPhi;
    // 6-4-E4. EZ = [EZ(x1), EZ(x2), ... , EZ(xn)]
    Ez = beta*X;
    // 6-4-E5. EZZsum = sum(EZZ(xi))
    EzzSum = n*(eyeK-beta*Lold) + beta*X*X.t()*beta.t();

    // 6-4-M1. update Lambda : Lnew
    Lnew = (X*Ez.t())*pinv(EzzSum);
    // 6-4-M2. update Phi    : Pnew from Ptmp
    Ptmp = (X - Lnew*Ez)*X.t()/n;
    Pnew = Ptmp.diag();

    // 6-4-Update
    inctol = norm(Lnew-Lold,"fro");
    Lold = Lnew;
    Pold = Pnew;
    if (inctol < tolerance){
      break;
    }
  }

  // 6-5. MLE solution for Z
  Pinv = diagmat(1/Pold);
  arma::mat Z = arma::solve(Lold.t()*Pinv*Lold,Lold.t()*sqrt(Pinv)*X);

  // 6-6. return results
  return Rcpp::List::create(Rcpp::Named("L")=Lold,
                            Rcpp::Named("Z")=Z,
                            Rcpp::Named("Pvec")=Pold);
}
arma::mat v2aux_adjproj(arma::mat& X){
  arma::mat Q, R;
  arma::qr_econ(Q,R,X);
  return(Q);
}
