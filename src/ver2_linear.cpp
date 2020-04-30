#include <RcppArmadillo.h>
#include "ver2_computation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/*
 * Class Definition =========================================================
 *    ClLLproc
 *
 * Main Functions ===========================================================
 *    (01) dt_pca   : PCA
 *    (02) dt_fa    : FA
 *    (03) dt_spca  : Sparse PCA
 *    (04) dt_mds   : MDS
 *    (05) dt_lasso : LASSO
 *    (06) dt_enet  : Elastic Net
 *    (07) dt_lmds  : Landmark MDS
 */

// Auxiliary Functions ======================================================
class ClLLproc{
private:
  std::string mytype;
public:
  // ClLLproc : constructor
  ClLLproc(std::string type){
    this->mytype = type;
  };

  // ClLLproc : functions
  arma::mat MainFunc(const arma::mat& X){
    // prepare : parameter
    int N = X.n_rows;
    int P = X.n_cols;

    // prepare : operators
    arma::rowvec mymean(P);
    arma::mat    mymult(P,P);

    std::string su = mytype;
    std::transform(su.begin(), su.end(), su.begin(), ::tolower);

    if (su=="null"){           // case 1. null
      mymean.fill(0.0);
      mymult.eye();
    } else if (su=="center"){  // case 2. center
      mymean.fill(0.0);
      mymean = arma::mean(X, 0);
      mymult.eye();
    } else if (su=="scale"){   // case 3. scale
      mymean.fill(0.0);
      mymult.fill(0.0);
      for (int p=0;p<P;p++){
        mymult(p,p) = 1.0/arma::stddev(X.col(p));
      }
    } else if (su=="cscale"){   // case 4. cscale
      mymean = arma::mean(X, 0);
      mymult.fill(0.0);
      // 4-2. scale (1/sd) by column
      for (int p=0;p<P;p++){
        mymult(p,p) = 1.0/arma::stddev(X.col(p));
      }
    } else {
      // others-1. mean should be 'center'
      mymean = arma::mean(X, 0);
      if (su=="decorrelate"){   // case 5. decorrelate
        arma::mat covX = arma::cov(X);
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, covX);

        mymult = eigvec;
      } else if (su=="whiten"){ // case 6. whiten
        arma::mat covX = arma::cov(X);
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, covX);

        arma::mat weight = arma::diagmat(1.0/arma::sqrt(eigval));
        mymult = eigvec*weight;
      } else { // error
        Rcpp::stop("Inadmissible Error.");
      }
    }

    // wrap the results
    return(arma::join_cols(mymean, mymult)); // [mu; multiplier]
  };
  std::string GetType(){
    std::string su = mytype;
    std::transform(su.begin(), su.end(), su.begin(), ::tolower);
    return(su);
  };
};



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
  ClLLproc init(ptype);
  arma::mat premat    = init.MainFunc(X);
  arma::rowvec mymean = premat.row(0);
  arma::mat    mymult = premat.rows(1,P);
  std::string  mytype = init.GetType();

  // data prep
  arma:mat Xtmp;
  if (mytype=="null"){
    Xtmp = X;
  } else {
    Xtmp.set_size(N,P);
    for (int n=0;n<N;n++){
      Xtmp.row(n) = X.row(n) - mymean;
    }
    Xtmp = Xtmp*mymult;
  }


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
  ClLLproc init(ptype);
  arma::mat premat    = init.MainFunc(X);
  arma::rowvec mymean = premat.row(0);
  arma::mat    mymult = premat.rows(1,P);
  std::string  mytype = init.GetType();

  // data prep
  arma::mat pX;
  if (mytype=="null"){
    pX = X;
  } else {
    pX.set_size(N,P);
    for (int n=0;n<N;n++){
      pX.row(n) = X.row(n) - mymean;
    }
    pX  = pX*mymult;
  }
  arma::mat tpX = arma::trans(pX);

  // main computation ---------------------------------------------------------
  // 1. pass onto v2aux_fa
  Rcpp::List output = v2aux_fa(tpX,ndim,maxiter,tolerance);

  // 2. after works
  arma::mat Y = output["Z"];

  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("type")=mytype,
    Rcpp::Named("mean")=mymean,
    Rcpp::Named("multiplier")=mymult,
    Rcpp::Named("algtype")="linear"
  );

  arma::mat LHS  = tpX*pX;
  arma::mat RHS  = tpX*Y.t();
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


// 03. SPCA ===============================================================
// auxiliary functions for dt_spca
arma::vec dt_spca_rk1vec(arma::mat& X){
  int p = X.n_rows;
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, X);

  arma::vec output = eigvec.tail_cols(1);
  return(output);
}
arma::mat dt_spca_deflation(arma::mat& Sig, arma::vec& Vec){
  int p = Vec.n_elem;
  arma::mat term1 = Sig*(Vec*Vec.t())*Sig; // column-vector convention
  double term2 = arma::dot((Sig*Vec), Vec);
  arma::mat output = Sig - (term1/term2);
  return(output);
}
// main function for spca
// [[Rcpp::export]]
Rcpp::List dt_spca(const arma::mat& X, int ndim, std::string ptype, double mu, double rho, const double abstol, const double reltol, const int maxiter){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
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


  // preprocessing
  ClLLproc init(ptype);
  arma::mat premat    = init.MainFunc(X);
  arma::rowvec mymean = premat.row(0);
  arma::mat    mymult = premat.rows(1,P);
  std::string  mytype = init.GetType();

  // data prep
  arma::mat pX;
  if (mytype=="null"){
    pX = X;
  } else {
    pX.set_size(N,P);
    for (int n=0;n<N;n++){
      pX.row(n) = X.row(n) - mymean;
    }
    pX  = pX*mymult;
  }

  // main computation ---------------------------------------------------------
  // 1. prepare for name changes in ADMM's SPCA function.
  int numpc = ndim;
  arma::mat Sigma = arma::cov(pX);
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

  // 3.   compute others
  // 3-1. trfinfo
  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("type")=mytype,
    Rcpp::Named("mean")=mymean,
    Rcpp::Named("multiplier")=mymult,
    Rcpp::Named("algtype")="linear"
  );
  // 3-2. update basis and Y
  arma::mat proj = v2aux_adjproj(basis);
  arma::mat Y    = pX*proj;



  // wrap and report ---------------------------------------------------
  // trfinfo
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("trfinfo") = info,
      Rcpp::Named("projection") = proj
  ));
}

// 04. MDS ================================================================
// [[Rcpp::export]]
Rcpp::List dt_mds(const arma::mat& X, int ndim, std::string ptype){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.mds : 'ndim' should be in [1,ncol(X)).");
  }

  // preprocessing
  ClLLproc init(ptype);
  arma::mat premat    = init.MainFunc(X);
  arma::rowvec mymean = premat.row(0);
  arma::mat    mymult = premat.rows(1,P);
  std::string  mytype = init.GetType();

  // data prep
  arma::mat pX;
  if (mytype=="null"){
    pX = X;
  } else {
    pX = X.each_row() - mymean;
    pX = pX*mymult;
  }
  arma::mat tpX = arma::trans(pX);

  // trfinfo
  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("type")=mytype,
    Rcpp::Named("mean")=mymean,
    Rcpp::Named("multiplier")=mymult,
    Rcpp::Named("algtype")="linear"
  );

  // main computation ---------------------------------------------------------
  // 1. eigendecomposition
  arma::mat U, V; arma::vec s;
  arma::svd(U,s,V,tpX);
  // 2. compute Y first
  arma::mat Y = V.head_cols(ndim)*arma::diagmat(s.head(ndim));
  // 3. compute Projection
  arma::mat LHS  = tpX*pX;
  arma::mat RHS  = tpX*Y;
  arma::mat proj = v2aux_solproj(LHS,RHS);

  // wrap and report ---------------------------------------------------
  // trfinfo
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("trfinfo") = info,
      Rcpp::Named("projection") = proj
  ));
}


// 05. LASSO ==============================================================
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
Rcpp::List dt_lasso(const arma::mat& X, int ndim, std::string ptype,
                    const arma::vec& y, bool ycenter, double lambda){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.lasso : 'ndim' should be in [1,ncol(X)).");
  }
  if (lambda <= arma::datum::eps){
    throw std::invalid_argument("* do.lasso : 'lambda' should be a nonnegative real number.");
  }

  // preprocessing
  ClLLproc init(ptype);
  arma::mat premat    = init.MainFunc(X);
  arma::rowvec mymean = premat.row(0);
  arma::mat    mymult = premat.rows(1,P);
  std::string  mytype = init.GetType();

  // data prep
  arma::mat pX;
  if (mytype=="null"){
    pX = X;
  } else {
    pX = X.each_row() - mymean;
    pX = pX*mymult;
  }

  arma::vec py = y;
  if (ycenter==true){
    double yc = arma::mean(y);
    py = y - yc;
  }

  // trfinfo
  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("type")=mytype,
    Rcpp::Named("mean")=mymean,
    Rcpp::Named("multiplier")=mymult,
    Rcpp::Named("algtype")="linear"
  );

  // main computation ---------------------------------------------------------
  // 1. main run of LASSO
  arma::vec x      = admm_lasso(pX, py, lambda);
  arma::vec lscore = arma::abs(x);

  // 2. find index of largest elements in magnitude
  arma::uvec idxtmp = arma::sort_index(lscore, "descend");
  arma::uvec idxvec = idxtmp.head(ndim);
  arma::uvec idRvec = idxvec + 1;

  // 3. compute projection by indicator
  arma::mat proj = v2aux_fid2proj(P,ndim,idxvec);
  arma::mat Y    = pX*proj;

  // wrap and report ---------------------------------------------------
  // trfinfo
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("featidx") = idRvec,
      Rcpp::Named("trfinfo") = info,
      Rcpp::Named("projection") = proj
  ));
}


// 06. Elastic Net ========================================================
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
  const int m = A.n_rows;
  const int n = A.n_cols;
  arma::mat U;
  arma::vec onesN(n,fill::ones);
  U = arma::chol(A.t()*A + rho*arma::diagmat(onesN));
  return(U);
}
arma::vec admm_enet(arma::mat& A, arma::vec& b,  double lambda, double alpha, double reltol, double abstol, int maxiter, double rho){
  // 1. get parameters
  const int m = A.n_rows;
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
Rcpp::List dt_enet(const arma::mat& X, int ndim, std::string ptype,
                   const arma::vec& y, bool ycenter, double lambda1, double lambda2){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
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

  // preprocessing
  ClLLproc init(ptype);
  arma::mat premat    = init.MainFunc(X);
  arma::rowvec mymean = premat.row(0);
  arma::mat    mymult = premat.rows(1,P);
  std::string  mytype = init.GetType();

  // data prep
  arma::mat pX;
  if (mytype=="null"){
    pX = X;
  } else {
    pX = X.each_row() - mymean;
    pX = pX*mymult;
  }

  arma::vec py = y;
  if (ycenter==true){
    double yc = arma::mean(y);
    py = y - yc;
  }

  // trfinfo
  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("type")=mytype,
    Rcpp::Named("mean")=mymean,
    Rcpp::Named("multiplier")=mymult,
    Rcpp::Named("algtype")="linear"
  );

  // main computation ---------------------------------------------------------
  // 0. extra parameters in enet function
  double lambda = 2*lambda2 + lambda1;
  double alpha  = lambda1/lambda;
  double abstol = 1e-4;
  double reltol = 1e-2;
  double rho = 1.0;
  int maxiter= 1000;


  // 1. main run of ENET
  arma::vec x      = admm_enet(pX, py, lambda, alpha, reltol, abstol, maxiter, rho);
  arma::vec lscore = arma::abs(x);

  // 2. find index of largest elements in magnitude
  arma::uvec idxtmp = arma::sort_index(lscore, "descend");
  arma::uvec idxvec = idxtmp.head(ndim);
  arma::uvec idRvec = idxvec + 1;

  // 3. compute projection by indicator
  arma::mat proj = v2aux_fid2proj(P,ndim,idxvec);
  arma::mat Y    = pX*proj;

  // wrap and report ---------------------------------------------------
  // trfinfo
  return(Rcpp::List::create(
      Rcpp::Named("Y") = Y,
      Rcpp::Named("featidx") = idRvec,
      Rcpp::Named("trfinfo") = info,
      Rcpp::Named("projection") = proj
  ));
}





// 07. LMDS ===================================================================
// [[Rcpp::export]]
Rcpp::List dt_lmds(const arma::mat& X, int ndim, std::string ptype, int npts){
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

  // preprocessing
  ClLLproc init(ptype);
  arma::mat premat    = init.MainFunc(X);
  arma::rowvec mymean = premat.row(0);
  arma::mat    mymult = premat.rows(1,P);
  std::string  mytype = init.GetType();

  // data prep
  arma::mat pX;
  if (mytype=="null"){
    pX = X;
  } else {
    pX.set_size(N,P);
    for (int n=0;n<N;n++){
      pX.row(n) = X.row(n) - mymean;
    }
    pX = pX*mymult;
  }


  // main computation ---------------------------------------------------------
  // 1. random selection
  arma::uvec idselect = arma::randperm(N,npts);
  // 2. do the pca
  arma::mat subX = pX.rows(idselect);
  arma::mat psdX = arma::cov(subX);
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, psdX); // ascending order issue
  arma::vec vars = arma::reverse(eigval.tail(ndim));
  arma::mat proj = arma::fliplr(eigvec.tail_cols(ndim));

  // 3. computation
  arma::mat Y = pX*proj;

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
      Rcpp::Named("projection") = proj,
      Rcpp::Named("trfinfo") = info
  ));
}
