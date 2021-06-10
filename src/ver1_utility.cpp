#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// (01) oos_linproj




// (01) oos_linproj  ===========================================================
// [[Rcpp::export]]
arma::mat oos_linproj(arma::mat& Xold, arma::mat& Yold, arma::mat& Xnew){
  // parameters
  // int n = Xold.n_rows;
  // int p = Xold.n_cols;
  // int m = Xnew.n_rows;
  // int ndim = Yold.n_cols;

  // // column means
  // arma::rowvec Xold_mean = arma::mean(Xold, 0);
  // arma::rowvec Yold_mean = arma::mean(Yold, 0);
  //
  // // subtract means
  // arma::mat Xtmp(n,p,fill::zeros);
  // arma::mat Ytmp(n,ndim,fill::zeros);
  // for (int i=0; i<n; i++){
  //   Xtmp.row(i) = Xold.row(i)-Xold_mean;
  //   Ytmp.row(i) = Yold.row(i)-Yold_mean;
  // }

  // solve pseudo
  arma::mat LHS = arma::trans(Xold)*Xold;
  arma::mat RHS = arma::trans(Xold)*Yold;
  arma::mat soltmp = arma::pinv(LHS)*RHS;

  // QR decomposition
  arma::mat Q;
  arma::mat R;
  arma::qr_econ(Q, R, soltmp);

  // output
  arma::mat Ynew = Xnew*Q;
  return(Ynew);
}
