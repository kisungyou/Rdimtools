#include <RcppArmadillo.h>
#include "ver1_computation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat v2aux_pca(arma::mat& X, int ndim){
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
arma::mat v2aux_pdist(arma::mat& X){
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
Rcpp::List v2aux_fa(arma::mat& X, const int k, const int maxiter, const double tolerance){
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
arma::mat v2aux_solproj(arma::mat& LHS, arma::mat& RHS){
  arma::mat Q, R;
  arma::mat X = arma::solve(LHS,RHS);
  arma::qr_econ(Q,R,X);
  return(Q);
}


arma::vec spca_gamma(arma::vec sigma, double r){
  const int p = sigma.n_elem;
  int indj = 0;
  double term1 = 0.0;
  double term2 = 0.0;
  for (int j=0;j<p;j++){
    term1 = sigma(j);
    for (int k=j;k<p;k++){
      term2 += sigma(k);
    }
    term2 = (term2-r)/(p-j);
    if (term1 > term2){
      indj = j;
      break;
    }
  }
  double theta = 0.0;
  for (int j=indj;j<p;j++){
    theta += sigma(j);
  }
  theta = (theta-r)/(p-indj);

  arma::vec output(p,fill::zeros);
  for (int i=0;i<p;i++){
    term1 = sigma(i)-theta;
    if (term1>0){
      output(i) = term1;
    }
  }
  return(output);
}

arma::mat spca_shrinkage(arma::mat A, const double tau){
  const int n = A.n_rows;
  arma::mat output(n,n,fill::zeros);
  double zij    = 0.0;
  double abszij = 0.0;
  double signer = 0.0;
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      zij = A(i,j);
      if (zij >= 0){
        signer = 1.0;
        abszij = zij;
      } else {
        signer = -1.0;
        abszij = -zij;
      }

      if (abszij > tau){
        output(i,j) = signer*(abszij-tau);
      }
    }
  }
  return(output);
}
arma::mat v2aux_spca(arma::mat& Sigma, const double reltol, const double abstol, const int maxiter, double mu, double rho){
  // 1. get parameters
  int p = Sigma.n_cols;

  // 2. set updating objects
  arma::mat Xold(p,p,fill::zeros);
  arma::mat Xnew(p,p,fill::zeros);
  arma::mat Yold(p,p,fill::zeros);
  arma::mat Ynew(p,p,fill::zeros);
  arma::mat Lold(p,p,fill::zeros);
  arma::mat Lnew(p,p,fill::zeros);

  arma::mat costX(p,p,fill::zeros);
  arma::mat costY(p,p,fill::zeros);
  arma::vec eigval(p,fill::zeros);   // for EVD of costX
  arma::mat eigvec(p,p,fill::zeros);

  // 3. iteration records
  arma::vec h_r_norm(maxiter,fill::zeros);
  arma::vec h_s_norm(maxiter,fill::zeros);
  arma::vec h_eps_pri(maxiter,fill::zeros);
  arma::vec h_eps_dual(maxiter,fill::zeros);

  // 4. main iteration
  int k=0;
  double ythr = mu*rho;
  double normX = 0.0;
  double normY = 0.0;
  for (k=0;k<maxiter;k++){
    // 4-1. update 'X'
    costX = Yold + mu*Lold + mu*Sigma;
    eig_sym(eigval, eigvec, costX);
    arma::vec gamma = spca_gamma(eigval, 1.0);
    Xnew = eigvec*arma::diagmat(gamma)*eigvec.t();

    // 4-2. update 'Y'
    costY = Xnew-mu*Lold;
    Ynew  = spca_shrinkage(costY, ythr);

    // 4-3. update 'L'
    Lnew  = Lold - (Xnew-Ynew)/mu;

    // 4-4. diagnostics for reporting
    h_r_norm(k) = arma::norm(Xnew-Ynew,"fro");
    h_s_norm(k) = arma::norm(Yold-Ynew,"fro")/mu;

    normX = arma::norm(Xnew,"fro");
    normY = arma::norm(Ynew,"fro");
    if (normX >= normY){
      h_eps_pri(k) = p*abstol + reltol*normX;
    } else {
      h_eps_pri(k) = p*abstol + reltol*normY;
    }
    h_eps_dual(k) = p*abstol + reltol*arma::norm(Lnew,"fro");

    // 4-5. updating and termination
    Xold = Xnew;
    Yold = Ynew;
    Lold = Lnew;

    if ((h_r_norm(k) < h_eps_pri(k))&&(h_s_norm(k)<h_eps_dual(k))){
      break;
    }
  }

  // 5. report results
  return(Xold); // ADMM's X object
}
arma::umat v2aux_knn(arma::mat& X, int k){
  int n = X.n_rows;
  int d = X.n_cols;
  arma::umat indices(n,k,fill::zeros);
  arma::vec  distvec(n,fill::zeros);
  arma::uvec sorted;
  arma::rowvec xi(d, fill::zeros);
  arma::rowvec xj(d, fill::zeros);

  for (int i=0;i<n;i++){
    distvec.fill(0.0);
    xi = X.row(i);
    for (int j=0;j<n;j++){
      xj = X.row(j);
      if (i==j){
        distvec(j) = 0.0;
      } else {
        distvec(j) = arma::norm(xi-xj,2);
      }
    }
    sorted = arma::sort_index(distvec,"ascend");
    indices.row(i) = arma::trans(sorted.rows(1,k));
  }
  return(indices);
}
arma::mat v2aux_fid2proj(int p, int ndim, arma::uvec idxvec){
  // replicate 'aux.featureindicator'
  // (p-by-ndim) indicator matrix for projection
  arma::mat output(p,ndim,fill::zeros);
  for (int i=0;i<ndim;i++){
    output(idxvec(i),i) = 1.0;
  }
  return(output);
}
// [[Rcpp::export]]
arma::vec v2aux_pagerank(arma::mat& A){
  // parameters
  int N = A.n_rows;
  double NN = static_cast<double>(N);
  double d  = 0.85; // standard choice

  int maxiter = 100;
  double thr  = 0.001;

  // prepare
  arma::vec sumA = arma::sum(A,1);
  arma::vec Kvec(N,fill::zeros);
  for (int n=0;n<N;n++){
    if (sumA(n) > 0){
      Kvec(n) = 1.0/sumA(n);
    }
  }
  arma::mat Kinv = arma::diagmat(Kvec);
  arma::mat M    = arma::trans(Kinv*A);

  // iterate
  arma::vec scterm(N,fill::zeros);
  arma::vec Rold(N,fill::zeros);
  arma::vec Rnew(N,fill::zeros);
  for (int n=0;n<N;n++){
    Rold(n)   = 1.0/NN;
    scterm(n) = ((1.0-d)/NN);
  }
  double Rinc = 0.0;
  for (int it=0;it<maxiter;it++){
    // update
    Rnew = d*M*Rold + scterm;
    Rinc = arma::norm(Rnew-Rold,2);
    Rold = Rnew;
    if ((Rinc < thr)&&(it>5)){
      break;
    }
  }

  // return
  return(Rold);
}

// [[Rcpp::export]]
arma::mat v2aux_pdist2(arma::mat &X, arma::mat &Y){
  int M = X.n_rows;
  int N = Y.n_rows;

  arma::mat output(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      output(m,n) = arma::norm(X.row(m)-Y.row(n),2);
    }
  }
  return(output);
}
