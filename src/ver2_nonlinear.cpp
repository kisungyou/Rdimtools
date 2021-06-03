#include <RcppArmadillo.h>
#include "ver1_computation.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/*
* Class Definition =========================================================
*    ClNLproc
*
* Main Functions ===========================================================
*    (01) dt_rpca   : RPCA
*    (03) dt_phate  : PHATE by Smita Krishnaswamy
*    (04) dt_mmds   : metric MDS
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




// (03) dt_phate ===============================================================
double dt_phate_entropy(arma::vec eigvals, int p){
  double threshold  = arma::datum::eps*100;
  arma::vec eigp    = arma::pow(eigvals, static_cast<double>(p));
  arma::vec normvec = eigp/arma::accu(eigp);

  double output = 0.0;
  double tgt    = 0.0;
  for (int i=0; i<eigvals.n_elem; i++){
    tgt = normvec(i);
    if (tgt > threshold){
      output += -tgt*std::log(tgt);
    }
  }
  return(output);
}
double dt_phate_regression(arma::vec x, arma::vec y){
  int N = x.n_elem;
  double xbar = arma::as_scalar(arma::mean(x));
  double ybar = arma::as_scalar(arma::mean(y));

  double beta1 = arma::dot(x-xbar, y-ybar);
  double beta2 = arma::accu(arma::pow(x-xbar, 2.0));

  double hat_beta  = beta1/beta2;
  double hat_alpha = ybar - hat_beta*xbar;

  double tgt = 0.0;
  double mse = 0.0;
  for (int n=0; n<N; n++){
    tgt  = (y(n) - (hat_alpha + hat_beta*x(n)));
    mse += tgt*tgt;
  }
  return(mse);
}
int dt_phate_linfits(arma::vec entropy){
  int N = entropy.n_elem;
  arma::vec mse(N,fill::zeros);
  mse(0)   = arma::datum::inf;
  mse(N-1) = arma::datum::inf;

  arma::vec vecx = arma::regspace<arma::vec>(1, 1, N);

  arma::vec vecx1(2,fill::zeros);
  arma::vec vecy1(2,fill::zeros);
  arma::vec vecx2(2,fill::zeros);
  arma::vec vecy2(2,fill::zeros);

  for (int i=1; i<(N-1); i++){
    // first fit
    vecx1.reset(); vecx1 = vecx.head(i+1);
    vecy1.reset(); vecy1 = entropy.head(i+1);

    // second fit
    vecx2.reset(); vecx2 = vecx.tail(N-(i+1));
    vecy2.reset(); vecy2 = entropy.tail(N-(i+1));

    // sum of two mse
    mse(i) = dt_phate_regression(vecx1,vecy1) + dt_phate_regression(vecx2,vecy2);
  }
  int idmin = mse.index_min();
  return(idmin + 1); // return the times
}
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
Rcpp::List dt_phate(arma::mat& X, int ndim, std::string ptype, int k, double alpha, std::string dtype, int maxiter, double abstol, bool smacof){
  // preliminary --------------------------------------------------------------
  // parameters
  int N = X.n_rows;
  if ((ndim < 1)||(ndim >= X.n_cols)){
    throw std::invalid_argument("* do.phate : 'ndim' should be in [1,ncol(X)).");
  }
  // preprocessing
  // ClNLproc init(ptype);
  // arma::mat premat    = init.MainFunc(X);
  // arma::rowvec mymean = premat.row(0);
  // arma::mat    mymult = premat.rows(1,X.n_cols);
  // std::string  mytype = init.GetType();
  arma::rowvec mymean(X.n_cols, fill::zeros);
  arma::mat mymult(X.n_cols, X.n_cols, fill::eye);
  std::string mytype = "null";

  // computation ---------------------------------------------------------------
  int i, j;
  // 1. pairwise distance computation : these are modified by nbdist
  arma::mat D(N,N,fill::zeros);
  for (i=0; i<(N-1); i++){
    for (j=(i+1); j<N; j++){
      D(i,j) = arma::norm(X.row(i)-X.row(j),2);
      D(j,i) = D(i,j);
    }
  }
  // 2. nearest distances
  arma::vec tgt(N,fill::zeros);
  arma::vec ND(N,fill::zeros);
  for (i=0; i<N; i++){
    tgt   = arma::sort(D.col(i),"ascend");
    ND(i) = tgt(k+1);
  }
  // 3. compute kernel matrix, rowsum, and markov transition matrix
  arma::mat K(N,N,fill::ones);
  for (i=0; i<(N-1); i++){
    for (j=(i+1); j<N; j++){
      K(i,j) = 0.5*(std::exp(-std::pow((D(i,j)/ND(i)), alpha)) + std::exp(-std::pow((D(i,j)/ND(j)), alpha)));
      K(j,i) = K(i,j);
    }
  }
  arma::vec Ksum = arma::sum(K, 1);
  arma::mat P = arma::diagmat(1.0/Ksum)*K;
  arma::mat A = arma::diagmat(1.0/arma::sqrt(Ksum))*K*arma::diagmat(1.0/arma::sqrt(Ksum));

  // 4. compute von-neumann entropy for best time-marching step
  arma::vec eigvals = arma::eig_sym(A);
  arma::vec  vec_entropy(200, fill::zeros);
  for (i=0; i<200; i++){
    vec_entropy(i) = dt_phate_entropy(eigvals, (i+1));
  }
  int times = dt_phate_linfits(vec_entropy);

  // 5. time marching
  arma::mat PT = P;
  if (times > 1){
    for (int tt=0; tt<(times-1); tt++){
      PT = PT*P;
    }
  }
  for (i=0; i<N; i++){
    PT.row(i) /= arma::accu(PT.row(i));
  }
  // 6. compute potential distance
  arma::mat PotDist(N,N,fill::zeros);
  arma::rowvec PTi(N,fill::zeros);
  arma::rowvec PTj(N,fill::zeros);
  for (i=0; i<(N-1); i++){
    PTi = PT.row(i);
    for (j=(i+1); j<N; j++){
      PTj = PT.row(j);
      if (dtype=="log"){
        PotDist(i,j) = arma::norm(arma::log(PTi)-arma::log(PTj), 2);
      } else if (dtype=="sqrt"){ // Riemannian failed
        PotDist(i,j) = arma::norm(arma::sqrt(PTi)-arma::sqrt(PTj),2);
      } else {
        PotDist(i,j) = arma::norm(PTi-PTj,2);
      }
      PotDist(j,i) = PotDist(i,j);
    }
  }
  PotDist *= arma::accu(D)/(static_cast<double>(N*N)*PotDist.max()); // arbitrary normalization for numerical reason

  // 7. apply MDS - either Classical or Metric
  arma::mat embed(N,ndim,fill::zeros);
  if (smacof==true){
    embed = dt_phate_mmds(PotDist, ndim, maxiter, abstol);
  } else {
    embed = dt_phate_cmds(PotDist, ndim);
  }


  // wrap and report ---------------------------------------------------
  // trfinfo
  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("type")=mytype,
    Rcpp::Named("mean")=mymean,
    Rcpp::Named("multiplier")=mymult,
    Rcpp::Named("algtype")="nonlinear"
  );
  return(Rcpp::List::create(
      Rcpp::Named("Y") = embed,
      Rcpp::Named("trfinfo") = info
  ));
}

