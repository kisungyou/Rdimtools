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
 *    dt_fa  : FA
 *    dt_pca : PCA
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
  ClLLproc init(ptype);
  arma::mat premat    = init.MainFunc(X);
  arma::rowvec mymean = premat.row(0);
  arma::mat    mymult = premat.rows(1,P);
  std::string  mytype = init.GetType();

  // data prep
  arma::mat pX(N,P,fill::zeros);
  arma::mat tpX(P,N,fill::zeros);
  for (int n=0;n<N;n++){
    pX.row(n) = X.row(n) - mymean;
  }
  pX  = pX*mymult;
  tpX = arma::trans(pX);

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
