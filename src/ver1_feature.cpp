#include <RcppArmadillo.h>
#include "ver1_computation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


// (uvec) dt_feature_smallidx   : return the k smallest indices
// (mat)  dt_feature_projection : best index to the projection matrix

// 01. CSCORE
// 02. LASSO
// 03. ENET

arma::uvec dt_feature_smallidx(arma::vec x, int k){
  arma::uvec indices = arma::sort_index(x, "ascend");
  arma::uvec topkvec = indices.head(k);
  return(topkvec);
}
arma::mat dt_feature_projection(int p, int ndim, arma::uvec idx){
  arma::mat output(p,ndim,fill::zeros);
  for (int i=0; i<ndim; i++){
    output(idx(i),i) = 1.0;
  }
  return(output);
}

// 01. CSCORE ==================================================================
arma::vec dt_cscore_scoresum(arma::mat X, arma::mat S){
  // parameter
  int N = X.n_rows;
  int P = X.n_cols;

  double tmpval = 0.0;
  double fdiff  = 0.0;
  arma::vec output(P,fill::zeros);
  for (int r=0;r<P;r++){
    tmpval = 0.0;
    for (int i=0;i<(N-1);i++){
      for (int j=(i+1);j<N;j++){
        fdiff = X(i,r) - X(j,r);
        tmpval += 2.0*(fdiff*fdiff)*S(i,j);
      }
    }
    output(r) = tmpval;
  }
  return(output);
}
// [[Rcpp::export]]
Rcpp::List dt_cscore(arma::mat& X, int ndim, arma::uvec& label, std::string myscore, double mylbd){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.lmds : 'ndim' should be in [1,ncol(X)).");
  }

  // computation ---------------------------------------------------------------
  // matrix for denoting same or different label
  arma::mat matSC(N,N,fill::zeros);
  arma::mat matSM(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      if (label(i)==label(j)){
        matSM(i,j) = 1;
        matSM(j,i) = 1;
      } else {
        matSC(i,j) = 1;
        matSC(j,i) = 1;
      }
    }
  }

  // compute elementary vectors
  arma::vec vecM = dt_cscore_scoresum(X, matSM);
  arma::vec vecC = dt_cscore_scoresum(X, matSC);

  // score according to the score type
  arma::vec rankvec;
  if (myscore=="ratio"){
    rankvec = vecM/vecC;
  } else {
    rankvec = vecM - mylbd*vecC;
  }

  // select the smallest ones
  arma::uvec idxvec     = dt_feature_smallidx(rankvec, ndim);
  arma::mat  projection = dt_feature_projection(P, ndim, idxvec);

  // wrap and report ---------------------------------------------------
  return(Rcpp::List::create(
      Rcpp::Named("Y") = X*projection,
      Rcpp::Named("cscore") = rankvec,
      Rcpp::Named("featidx") = idxvec+1,
      Rcpp::Named("projection") = projection,
      Rcpp::Named("algorithm") = "linear:CSCORE"
  ));
}

// 02. LASSO ===================================================================
arma::vec lasso_shrinkage(arma::vec a, const double kappa){
  const int n = a.n_elem;
  arma::vec y(n,fill::zeros);
  for (int i=0;i<n;i++){
    // first term : max(0, a-kappa)
    if (a(i)-kappa > 0){
      y(i) = a(i)-kappa;
    }
    // second term : -max(0, -a-kappa)
    if (-a(i)-kappa > 0){
      y(i) = y(i) + a(i) + kappa;
    }
  }
  return(y);
}
double lasso_objective(arma::mat A, arma::vec b, const double lambda, arma::vec x, arma::vec z){
  return(norm(A*x-b,2)/2 + lambda*norm(z,1));
}
arma::mat lasso_factor(arma::mat A, double rho){
  const int m = A.n_rows;
  const int n = A.n_cols;
  arma::mat U;
  if (m>=n){ // skinny case
    arma::vec onesN(n,fill::ones);
    U = chol(A.t()*A + rho*diagmat(onesN));
  } else {
    arma::vec onesM(m,fill::ones);
    U = chol(diagmat(onesM)+(1.0/rho)*(A*A.t()));
  }
  return(U);
}
// [[Rcpp::export]]
arma::vec admm_lasso(arma::mat& A, arma::vec& b, double lambda){
  // 0. default values
  double abstol = 1e-4;
  double reltol = 1e-2;
  double maxiter = 1000;
  double rho   = 1.0;
  double alpha = 1.0;

  // 1. get parameters
  const int m = A.n_rows;
  const int n = A.n_cols;

  // 2. set ready
  arma::vec xinit(n,fill::randu);
  arma::vec x(n,fill::zeros);
  arma::vec z(n,fill::zeros);
  arma::vec u(n,fill::zeros);
  arma::vec q(n,fill::zeros);
  arma::vec zold(n,fill::zeros);
  arma::vec x_hat(n,fill::zeros);

  // 3. precompute static variables for x-update and factorization
  arma::mat Atb = A.t()*b;
  arma::mat U   = lasso_factor(A,rho); // returns upper
  arma::mat L   = U.t();

  // 4. iteration
  arma::vec h_objval(maxiter,fill::zeros);
  arma::vec h_r_norm(maxiter,fill::zeros);
  arma::vec h_s_norm(maxiter,fill::zeros);
  arma::vec h_eps_pri(maxiter,fill::zeros);
  arma::vec h_eps_dual(maxiter,fill::zeros);

  double rho2  = rho*rho;
  double sqrtn = std::sqrt(static_cast<float>(n));
  int k;
  for (k=0;k<maxiter;k++){
    // 4-1. update 'x'
    q = Atb + rho*(z-u); // temporary value
    if (m>=n){
      x = solve(trimatu(U),solve(trimatl(L),q));
    } else {
      x = q/rho - (A.t()*solve(trimatu(U),solve(trimatl(L),A*q)))/rho2;
    }

    // 4-2. update 'z' with relaxation
    zold = z;
    x_hat = alpha*x + (1 - alpha)*zold;
    z = lasso_shrinkage(x_hat + u, lambda/rho);

    // 4-3. update 'u'
    u = u + (x_hat - z);

    // 4-3. dianostics, reporting
    h_objval(k) = lasso_objective(A,b,lambda,x,z);
    h_r_norm(k) = norm(x-z);
    h_s_norm(k) = norm(-rho*(z-zold));
    if (norm(x)>norm(-z)){
      h_eps_pri(k) = sqrtn*abstol + reltol*norm(x);
    } else {
      h_eps_pri(k) = sqrtn*abstol + reltol*norm(-z);
    }
    h_eps_dual(k) = sqrtn*abstol + reltol*norm(rho*u);

    // 4-4. termination
    if ((h_r_norm(k) < h_eps_pri(k))&&(h_s_norm(k)<h_eps_dual(k))){
      break;
    }
  }

  // 5. report results
  return(x);
}
// [[Rcpp::export]]
Rcpp::List dt_lasso(arma::mat& X, int ndim, arma::vec& y, double lambda){
  // preliminary --------------------------------------------------------------
  // parameters
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.lasso : 'ndim' should be in [1,ncol(X)).");
  }
  if (lambda <= arma::datum::eps){
    throw std::invalid_argument("* do.lasso : 'lambda' should be a nonnegative real number.");
  }

  // main computation ---------------------------------------------------------
  // 1. main run of LASSO
  arma::vec x      = admm_lasso(X, y, lambda);
  arma::vec lscore = arma::abs(x);

  // 2. find index of largest elements in magnitude
  arma::uvec idxtmp = arma::sort_index(lscore, "descend");
  arma::uvec idxvec = idxtmp.head(ndim);
  arma::uvec idRvec = idxvec + 1;

  // 3. compute projection by indicator
  arma::mat proj = v2aux_fid2proj(P,ndim,idxvec);
  arma::mat Y    = X*proj;

  // wrap and report ---------------------------------------------------
  // trfinfo
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("featidx") = idRvec,
      Rcpp::Named("projection") = proj,
      Rcpp::Named("algorithm") = "linear:LASSO"
  ));
}


// 03. Elastic Net ========================================================
arma::vec enet_shrinkage(arma::vec a, double kappa){
  const int n = a.n_elem;
  arma::vec y(n,fill::zeros);
  for (int i=0;i<n;i++){
    // first term : max(0, a-kappa)
    if (a(i)-kappa > 0){
      y(i) = a(i)-kappa;
    }
    // second term : -max(0, -a-kappa)
    if (-a(i)-kappa > 0){
      y(i) = y(i) + a(i) + kappa;
    }
  }
  return(y);
}
double enet_objective(arma::mat& A, arma::vec& b, double lambda, double alpha, arma::vec& x, arma::vec& z){
  //return(pow(norm(A*x-b,2),2)/2 + lambda*alpha*norm(z,1) + 0.5*(1-alpha)*lambda*pow(norm(x,2),2));
  return(norm(A*x-b,2)/2+lambda*alpha*norm(z,1)+0.5*(1-alpha)*lambda*norm(x,2));
}
arma::mat enet_factor(arma::mat& A, double rho){
  const int n = A.n_cols;
  arma::mat U;
  arma::vec onesN(n,fill::ones);
  U = arma::chol(A.t()*A + rho*arma::diagmat(onesN));
  return(U);
}
arma::vec admm_enet(arma::mat& A, arma::vec& b,  double lambda, double alpha, double reltol, double abstol, int maxiter, double rho){
  // 1. get parameters
  const int n = A.n_cols;
  double gamma = lambda*(1-alpha)+rho;

  // 2. set ready
  arma::vec x(n,fill::zeros);
  arma::vec z(n,fill::zeros);
  arma::vec u(n,fill::zeros);
  arma::vec q(n,fill::zeros);
  arma::vec zold(n,fill::zeros);
  arma::vec x_hat(n,fill::zeros);

  // 3. precompute static variables for x-update and factorization
  arma::mat Atb = A.t()*b;
  arma::mat U   = enet_factor(A,gamma); // returns upper
  arma::mat L   = U.t();


  // 4. iteration
  arma::vec h_objval(maxiter,fill::zeros);
  arma::vec h_r_norm(maxiter,fill::zeros);
  arma::vec h_s_norm(maxiter,fill::zeros);
  arma::vec h_eps_pri(maxiter,fill::zeros);
  arma::vec h_eps_dual(maxiter,fill::zeros);

  // double rho2 = rho*rho;

  double sqrtn = std::sqrt(static_cast<float>(n));
  int k;
  for (k=0; k<maxiter; k++){
    // 4-1. update 'x'
    q = Atb + rho*(z-u); // temporary value
    x = solve(trimatu(U),solve(trimatl(L),q));

    // 4-2. update 'z'
    zold = z;
    z = enet_shrinkage(x + u, lambda*alpha/rho);

    // 4-3. update 'u'
    u = u + x - z;

    // 4-3. dianostics, reporting
    h_objval(k) = enet_objective(A,b,lambda,alpha,x,z);
    h_r_norm(k) = arma::norm(x-z);
    h_s_norm(k) = arma::norm(-rho*(z-zold));
    if (norm(x)>norm(-z)){
      h_eps_pri(k) = sqrtn*abstol + reltol*norm(x);
    } else {
      h_eps_pri(k) = sqrtn*abstol + reltol*norm(-z);
    }
    h_eps_dual(k) = sqrtn*abstol + reltol*norm(rho*u);
    // 4-4. termination
    if ((h_r_norm(k) < h_eps_pri(k))&&(h_s_norm(k)<h_eps_dual(k))){
      break;
    }
  }
  // 5. report results
  return(x);
}

// [[Rcpp::export]]
Rcpp::List dt_enet(arma::mat& X, int ndim, arma::vec& y, double lambda1, double lambda2){
  // preliminary --------------------------------------------------------------
  // parameters
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.enet : 'ndim' should be in [1,ncol(X)).");
  }
  if (lambda1 <= arma::datum::eps){
    throw std::invalid_argument("* do.enet : 'lambda1' should be a nonnegative real number.");
  }
  if (lambda2 <= arma::datum::eps){
    throw std::invalid_argument("* do.enet : 'lambda2' should be a nonnegative real number.");
  }

  // main computation ---------------------------------------------------------
  // 0. extra parameters in enet function
  double lambda = 2*lambda2 + lambda1;
  double alpha  = lambda1/lambda;
  double abstol = 1e-4;
  double reltol = 1e-2;
  double rho = 1.0;
  int maxiter= 1000;


  // 1. main run of ENET
  arma::vec x      = admm_enet(X, y, lambda, alpha, reltol, abstol, maxiter, rho);
  arma::vec lscore = arma::abs(x);

  // 2. find index of largest elements in magnitude
  arma::uvec idxtmp = arma::sort_index(lscore, "descend");
  arma::uvec idxvec = idxtmp.head(ndim);
  arma::uvec idRvec = idxvec + 1;

  // 3. compute projection by indicator
  arma::mat proj = v2aux_fid2proj(P,ndim,idxvec);
  arma::mat Y    = X*proj;

  // wrap and report ---------------------------------------------------
  // trfinfo
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("featidx") = idRvec,
      Rcpp::Named("projection") = proj,
      Rcpp::Named("algorithm") = "linear:ENET"
  ));
}

