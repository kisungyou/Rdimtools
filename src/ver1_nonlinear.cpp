#include <RcppArmadillo.h>
#include "ver1_computation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// 01. MMDS
// 02. PHATE-partial
// 03. RPCA

// 01. MMDS : metric MDS -------------------------------------------------------
arma::mat dt_mmds_core(arma::mat D, int ndim, int maxiter, double abstol){
  // initialization with CMDS
  int N = D.n_rows;
  arma::mat D2 = arma::pow(D, 2.0);
  arma::mat J  = arma::eye<arma::mat>(N,N) - (arma::ones<arma::mat>(N,N)/(static_cast<double>(N)));
  arma::mat B  = -0.5*J*D2*J;  arma::vec eigval;  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, B);

  arma::mat old_X = eigvec.tail_cols(ndim)*arma::diagmat(arma::sqrt(eigval.tail(ndim)));
  arma::mat old_D(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      old_D(i,j) = arma::norm(old_X.row(i)-old_X.row(j),2);
      old_D(j,i) = old_D(i,j);
    }
  }

  arma::mat new_X(N,ndim,fill::zeros);
  arma::mat new_D(N,N,fill::zeros);
  double old_cost = 0.0;
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      old_cost += std::pow(D(i,j)-old_D(i,j), 2.0);
    }
  }
  double new_cost = 0.0;

  arma::mat BZ(N,N,fill::zeros);
  double dijZ   = 0.0;
  double epsthr = 100*arma::datum::eps;
  double inctol = 0.0;
  double diagsum = 0.0;

  for (int it=0; it<maxiter; it++){
    // compute updater B(Z); (Borg p155); Z=old_X
    BZ.fill(0.0);
    // off-diagonal first
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        // dijZ = arma::norm(old_X.row(i)-old_X.row(j),2);
        dijZ = old_D(i,j);
        if (dijZ > epsthr){
          BZ(i,j) = -D(i,j)/dijZ;
          BZ(j,i) = BZ(i,j);
        }
      }
    }
    // diagonal part
    for (int i=0; i<N; i++){
      diagsum = -arma::accu(BZ.row(i));
      BZ(i,i) = diagsum;
    }
    // updater
    new_X    = (BZ*old_X)/(static_cast<double>(N));
    new_D.fill(0.0);
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        new_D(i,j) = arma::norm(new_X.row(i)-new_X.row(j),2);
        new_D(j,i) = new_D(i,j);
      }
    }
    new_cost = 0.0;
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        new_cost += std::pow(D(i,j)-new_D(i,j), 2.0);
      }
    }

    inctol   = old_cost - new_cost;
    old_X    = new_X;
    old_cost = new_cost;
    if (inctol < abstol){
      break;
    }
  }
  return(old_X);
}
// [[Rcpp::export]]
Rcpp::List dt_mmds(arma::mat& X, int ndim, int maxiter, double abstol){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= X.n_cols)){
    throw std::invalid_argument("* do.mmds : 'ndim' should be in [1,ncol(X)).");
  }

  // computation ---------------------------------------------------------------
  // distance computation
  arma::mat D(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      D(i,j) = arma::norm(X.row(i)-X.row(j),2);
      D(j,i) = D(i,j);
    }
  }
  // compute
  arma::mat embed = dt_mmds_core(D, ndim, maxiter, abstol);

  // wrap and report ---------------------------------------------------
  return(Rcpp::List::create(
      Rcpp::Named("Y") = embed,
      Rcpp::Named("algorithm") = "nonlinear:MMDS"
  ));
}

// 02. PHATE-partial ===========================================================
arma::mat dt_phate_cmds(arma::mat pdist, int ndim){
  int N = pdist.n_rows;
  arma::mat D2 = arma::pow(pdist, 2.0);
  arma::mat J  = arma::eye<arma::mat>(N,N) - (arma::ones<arma::mat>(N,N)/(static_cast<double>(N)));
  arma::mat B  = -0.5*J*D2*J;

  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, B);
  arma::mat Y = eigvec.tail_cols(ndim)*arma::diagmat(arma::sqrt(eigval.tail(ndim)));
  return(Y);
}
arma::mat dt_phate_mmds(arma::mat D, int ndim, int maxiter, double abstol){
  // initialization with CMDS
  int N = D.n_rows;
  // arma::mat D2 = arma::pow(D, 2.0);
  // arma::mat J  = arma::eye<arma::mat>(N,N) - (arma::ones<arma::mat>(N,N)/(static_cast<double>(N)));
  // arma::mat B  = -0.5*J*D2*J;  arma::vec eigval;  arma::mat eigvec;
  // arma::eig_sym(eigval, eigvec, B);

  // arma::mat old_X = eigvec.tail_cols(ndim)*arma::diagmat(arma::sqrt(eigval.tail(ndim)));
  arma::mat old_X = dt_phate_cmds(D, ndim);
  arma::mat old_D(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      old_D(i,j) = arma::norm(old_X.row(i)-old_X.row(j),2);
      old_D(j,i) = old_D(i,j);
    }
  }

  arma::mat new_X(N,ndim,fill::zeros);
  arma::mat new_D(N,N,fill::zeros);
  double old_cost = 0.0;
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      old_cost += std::pow(D(i,j)-old_D(i,j), 2.0);
    }
  }
  double new_cost = 0.0;

  arma::mat BZ(N,N,fill::zeros);
  double dijZ   = 0.0;
  double epsthr = 100*arma::datum::eps;
  double inctol = 0.0;
  double diagsum = 0.0;

  for (int it=0; it<maxiter; it++){
    // compute updater B(Z); (Borg p155); Z=old_X
    BZ.fill(0.0);
    // off-diagonal first
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        // dijZ = arma::norm(old_X.row(i)-old_X.row(j),2);
        dijZ = old_D(i,j);
        if (dijZ > epsthr){
          BZ(i,j) = -D(i,j)/dijZ;
          BZ(j,i) = BZ(i,j);
        }
      }
    }
    // diagonal part
    for (int i=0; i<N; i++){
      diagsum = -arma::accu(BZ.row(i));
      BZ(i,i) = diagsum;
    }
    // updater
    new_X    = (BZ*old_X)/(static_cast<double>(N));
    new_D.fill(0.0);
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        new_D(i,j) = arma::norm(new_X.row(i)-new_X.row(j),2);
        new_D(j,i) = new_D(i,j);
      }
    }
    new_cost = 0.0;
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        new_cost += std::pow(D(i,j)-new_D(i,j), 2.0);
      }
    }

    inctol   = old_cost - new_cost;
    old_X    = new_X;
    old_cost = new_cost;
    if (inctol < abstol){
      break;
    }
  }
  return(old_X);
}
// [[Rcpp::export]]
Rcpp::List dt_phate_partial(arma::mat& P, int ndim, std::string dtype, int maxiter, double abstol, bool smacof){
  // get some parameters
  int N = P.n_rows;

  // automatic error catch-up with 'dtype' : log may incur numerical errors a lot.
  std::string mydtype = dtype;
  if (dtype=="log"){
    if (P.min() < arma::datum::eps){
      mydtype = "sqrt";
      Rcpp::warning("* do.phate : numerical underflow triggered by 'log' distance. Automatically change to 'sqrt' distance.");
    }
  }

  // compute pairwise distance given a row-stochastic matrix
  arma::mat PotDist(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      if (mydtype=="log"){
        PotDist(i,j) = arma::norm(arma::log(P.row(i))-arma::log(P.row(j)), 2);
      } else if (mydtype=="sqrt"){ // Riemannian failed
        PotDist(i,j) = arma::norm(arma::sqrt(P.row(i))-arma::sqrt(P.row(j)),2);
      } else {
        PotDist(i,j) = arma::norm(P.row(i)-P.row(j),2);
      }
      PotDist(j,i) = PotDist(i,j);
    }
  }

  // 7. apply MDS - either Classical or Metric
  arma::mat embed(N,ndim,fill::zeros);
  if (smacof==true){
    embed = dt_phate_mmds(PotDist, ndim, maxiter, abstol);
  } else {
    embed = dt_phate_cmds(PotDist, ndim);
  }

  // wrap and report ---------------------------------------------------
  return(Rcpp::List::create(
      Rcpp::Named("Y") = embed,
      Rcpp::Named("algorithm") = "nonlinear:PHATE"
  ));
}

// 03. RPCA ====================================================================
arma::vec shrink_vec_rpca(arma::vec x, double tau){
  const int n = x.n_elem;
  arma::vec output(n,fill::zeros);
  double xij    = 0.0;
  double absxij = 0.0;
  double signer = 0.0;
  for (int i=0;i<n;i++){
    xij = x(i);
    if (xij >= 0){
      signer = 1.0;
      absxij = xij;
    } else {
      signer = 1.0;
      absxij = -xij;
    }
    if (absxij > tau){
      output(i) = signer*(absxij-tau);
    }
  }
  return(output);
}
arma::mat shrink_mat_rpca(arma::mat A, const double tau){
  const int n = A.n_rows;
  const int p = A.n_cols;
  arma::mat output(n,p,fill::zeros);
  double zij    = 0.0;
  double abszij = 0.0;
  double signer = 0.0;
  for (int i=0;i<n;i++){
    for (int j=0;j<p;j++){
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
arma::mat rpca_vectorpadding(arma::vec x, const int n, const int p){
  arma::mat output(n,p,fill::zeros);
  if (n<p){
    for (int i=0;i<n;i++){
      output(i,i) = x(i);
    }
  } else {
    for (int j=0;j<p;j++){
      output(j,j) = x(j);
    }
  }
  return(output);
}
Rcpp::List admm_rpca(arma::mat& M, const double tol, const int maxiter,
                     double mu, double lambda){
  // 1. get parameters
  const int n1 = M.n_rows;
  const int n2 = M.n_cols;

  double invmu = 1/mu;
  double lbdmu = lambda/mu;

  // 2. set updating objects
  arma::mat Lold(n1,n2,fill::zeros);
  arma::mat Lnew(n1,n2,fill::zeros);
  arma::mat Sold(n1,n2,fill::zeros);
  arma::mat Snew(n1,n2,fill::zeros);
  arma::mat Yold(n1,n2,fill::zeros);
  arma::mat Ynew(n1,n2,fill::zeros);

  arma::mat costL(n1,n2,fill::zeros);
  arma::mat costS(n1,n2,fill::zeros);
  arma::mat costY(n1,n2,fill::zeros);
  arma::mat spadding(n1,n2,fill::zeros);

  arma::mat svdU;
  arma::vec svds;
  arma::mat svdV;
  arma::vec vecshrinkage;

  // 3. iteration records
  arma::vec vectolerance(maxiter,fill::zeros);

  // 4. main iteration
  int k=0;
  double norm1 = 0.0;                 // error LHS term
  double norm2 = arma::norm(M,"fro"); // error RHS term
  double normratio = 0.0;
  for (k=0;k<maxiter;k++){
    //  4-1. update L
    costL = (M - Sold + Yold*invmu);                   // compute term to be decomposed
    svd(svdU, svds, svdV, costL);                      // svd decomposition
    vecshrinkage = shrink_vec_rpca(svds, invmu);       // do shrinkage on singular vector
    spadding = rpca_vectorpadding(vecshrinkage,n1,n2); // we need zero padding on diagmat one
    Lnew = svdU*spadding*svdV.t();                     // update L

    // 4-2. update S
    costS = (M-Lnew+Yold*invmu);                // compute term to be shrinked
    Snew  = shrink_mat_rpca(costS, lbdmu);        // update S

    // 4-3. update Y
    Ynew  = Yold + mu*(M-Lnew-Snew);

    // 4-4. compute error
    norm1     = arma::norm(M-Lnew-Snew,"fro");
    normratio = norm1/norm2;
    vectolerance(k) = normratio;

    // 4-5. updating and termination
    Lold = Lnew;
    Sold = Snew;
    Yold = Ynew;

    if (normratio < tol){
      break;
    }
  }

  // 5. report results
  List output;
  output["L"] = Lold;
  output["S"] = Sold;
  output["k"] = k;
  output["errors"] = vectolerance;
  return(output);
}
// [[Rcpp::export]]
Rcpp::List dt_rpca(arma::mat& X, double mu, double lambda, int maxiter, double abstol){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  if (mu < arma::datum::eps){
    throw std::invalid_argument("* do.rpca : 'mu' should be a positive real number.");
  }
  if (lambda <= arma::datum::eps){
    throw std::invalid_argument("* do.rpca : 'lambda' should be a nonnegative real number.");
  }

  // main computation ---------------------------------------------------------
  // 1. pass to ADMM
  Rcpp::List admmrun = admm_rpca(X, abstol, maxiter, mu, lambda);
  // 2. select L and S
  arma::mat L = admmrun["L"];
  arma::mat S = admmrun["S"];


  // wrap and report ---------------------------------------------------
  return(Rcpp::List::create(
      Rcpp::Named("L") = L,
      Rcpp::Named("S") = S,
      Rcpp::Named("algorithm") = "nonlinear:RPCA"
  ));
}


