#ifndef RDIMTOOLS_VER2_CLASSES_H
#define RDIMTOOLS_VER2_CLASSES_H

#include <RcppArmadillo.h>

// ClPreproc ------------------------------------------------------
// first row is 'mean' for row vector,
// all the other rows are for 'multiplier'
class ClPreproc{
private:
  std::string mytype;
public:
  // ClPreproc : constructor
  ClPreproc(std::string type);

  // ClPreproc : functions
  arma::mat MainFunc(const arma::mat& X);
  std::string GetType();
};

// ClPreproc : constructor
ClPreproc::ClPreproc(std::string type){
  this->mytype = type;
}

// ClPreproc : functions
std::string ClPreproc::GetType(){
  std::string su = mytype;
  std::transform(su.begin(), su.end(), su.begin(), ::tolower);
  return(su);
}
arma::mat ClPreproc::MainFunc(const arma::mat& X){
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
}

#endif
