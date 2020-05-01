#include <RcppArmadillo.h>
#include "ver2_computation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/*
* Class Definition =========================================================
*    ClNLproc
*
* Main Functions ===========================================================
*    (01) dt_rpca : RPCA
*/

// Auxiliary Functions ======================================================
class ClNLproc{
private:
  std::string mytype;
public:
  // ClNLproc : constructor
  ClNLproc(std::string type){
    this->mytype = type;
  };

  // ClNLproc : functions
  arma::mat MainFunc(arma::mat& X){
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

// 01. RPCA ==================================================================
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
Rcpp::List dt_rpca(arma::mat& X, int ndim, std::string ptype, double mu, double lambda){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  if ((ndim < 1)||(ndim >= P)){
    throw std::invalid_argument("* do.rpca : 'ndim' should be in [1,ncol(X)).");
  }
  if (mu < arma::datum::eps){
    throw std::invalid_argument("* do.rpca : 'mu' should be a positive real number.");
  }
  if (lambda <= arma::datum::eps){
    throw std::invalid_argument("* do.rpca : 'lambda' should be a nonnegative real number.");
  }


  // preprocessing
  ClNLproc init(ptype);
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
  // 1. pass to ADMM
  Rcpp::List admmrun = admm_rpca(pX, 1e-7, 1000, mu, lambda);
  // 2. select L and S
  arma::mat L = admmrun["L"];
  arma::mat S = admmrun["S"];


  // wrap and report ---------------------------------------------------
  // trfinfo
  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("type")=mytype,
    Rcpp::Named("mean")=mymean,
    Rcpp::Named("multiplier")=mymult,
    Rcpp::Named("algtype")="nonlinear"
  );
  return(Rcpp::List::create(
      Rcpp::Named("L") = L,
      Rcpp::Named("S") = S,
      Rcpp::Named("trfinfo") = info
  ));
}
