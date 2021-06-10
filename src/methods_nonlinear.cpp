#include <RcppDist.h>
#include "methods_nonlinear.h"

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// subroutines
//    1. Compute Q for Low-Dimensional Values
arma::mat computeQ(arma::mat& Y){
  int n = Y.n_cols;
  double cdenom = 0.0;
  arma::mat Q(n,n,fill::zeros);
  for (int i=0;i<n;i++){
    cdenom = 0;
    for (int k=0;k<n;k++){
      if (i==k){
        cdenom += 0;
      } else{
        cdenom += std::exp(static_cast<float>(-pow(arma::norm(Y.col(i)-Y.col(k)),2.0)));
      }
    }
    for (int j=0;j<n;j++){
      if (i==j){
        Q(i,j) = 0;
      } else {
        Q(i,j) = std::exp(static_cast<float>(-pow(norm(Y.col(i)-Y.col(j)),2.0)))/cdenom;
      }
    }
  }
  return(Q);
}

// 1. SNE : Stochastic Neighbor Embedding
// [[Rcpp::export]]
arma::mat method_sne(arma::mat& P, int ndim0, double eta0,
                     int maxiter0, double jitter0, double decay0,
                     double momentum0){
  // all parameters
  int ndim = ndim0;
  int maxiter = maxiter0;

  double eta = eta0;
  double jitter = jitter0;
  double decay  = decay0;
  double momentum = momentum0;

  // 1-1. Initialize
  int n = P.n_cols;
  arma::mat Y = arma::mat(ndim,n,fill::randn)*0.0001;
  arma::mat dC(ndim,n,fill::zeros);
  arma::mat y_incs(ndim,n,fill::zeros);
  arma::mat mask = (1e-10)*arma::ones<arma::mat>(n,n);

  // P = arma::max(mask,P);
  arma::mat myP(n,n,fill::zeros);
  myP = arma::max(mask, P);

  // 1-2. Main Iteration
  arma::mat Q(n,n,fill::zeros);
  arma::mat PQ(n,n,fill::zeros);
  for (int it=0;it<maxiter;it++){
    dC.fill(0.0);
    Q = computeQ(Y);
    Q = arma::max(mask,Q);
    for (int i=0;i<n;i++){
      for (int j=0;j<n;j++){
        if (j!=i){
          dC.col(i) += 2.0*(myP(j,i)+myP(i,j)-Q(i,j)-Q(j,i))*(Y.col(i)-Y.col(j));
        }
      }
    }
    y_incs = momentum*y_incs - eta*dC;
    Y += y_incs;
    mat JitMat = jitter*randn<mat>(ndim,n);
    Y += JitMat;
    vec MeanY = arma::mean(Y,1);
    for (int i=0;i<n;i++){
      Y.col(i) -= MeanY;
    }
    jitter *= decay;
  }
  return(Y);
}


// 2. Symmetric SNE : Stochastic Neighbor Embedding
// [[Rcpp::export]]
arma::mat method_snesym(arma::mat& P, int ndim0, double eta0,
                     int maxiter0, double jitter0, double decay0,
                     double momentum0){
  // all parameters
  int ndim = ndim0;
  int maxiter = maxiter0;

  double eta = eta0;
  double jitter = jitter0;
  double decay  = decay0;
  double momentum = momentum0;


  // (1) Initialize
  int n = P.n_cols;
  arma::mat Y = 0.0001*randn<arma::mat>(ndim,n);
  arma::mat dC(ndim,n,fill::zeros);
  arma::mat y_incs(ndim,n);
  arma::mat mask = (1e-10)*ones<mat>(n,n);
  // myP = arma::max(mask, myP);

  arma::mat myP = (P + P.t())/(2.0*static_cast<double>(n));
  myP = arma::max(mask, myP);

  // arma::uvec myPsmall = arma::find(myP < 1e-10);
  // myP(myPsmall) = 1e-10;
  // P = (P+P.t());
  // P /= (2*n);
  // P = arma::max(mask,P);

  // (2) Main Iteration
  arma::mat Q(n,n);
  arma::mat PQ(n,n);
  for (int it=0;it<maxiter;it++){
    dC.zeros();
    Q = computeQ(Y);
    Q = (Q+Q.t())/2;
    Q = arma::max(mask,Q);
    for (int i=0;i<n;i++){
      for (int j=0;j<n;j++){
        if (j!=i){
          dC.col(i) += 4.0*(myP(i,j)-Q(i,j))*(Y.col(i)-Y.col(j));
        }
      }
    }
    y_incs = momentum*y_incs - eta*dC;
    Y += y_incs;
    mat JitMat = jitter*randn<mat>(ndim,n);
    Y += JitMat;
    arma::vec MeanY = mean(Y,1);
    for (int i=0;i<n;i++){
      Y.col(i) -= MeanY;
    }
    jitter *= decay;
  }
  return(Y);
}

// 3. tSNE : t-Stochastic Neighbor Embedding
// [[Rcpp::export]]
arma::mat method_tsne(arma::mat& P, int ndim0, double eta0,
                        int maxiter0, double jitter0, double decay0,
                        double momentum0){
  // all parameters
  int ndim = ndim0;
  int maxiter = maxiter0;

  double eta = eta0;
  double jitter = jitter0;
  double decay  = decay0;
  double momentum = momentum0;

  // 3-1. Initialize
  int n = P.n_cols;
  arma::mat Y = 0.0001*randn<mat>(ndim,n);
  arma::mat dC(ndim,n);
  arma::mat y_incs(ndim,n);
  arma::mat mask = (1e-10)*ones<mat>(n,n);

  arma::mat myP = (P + P.t())/(2.0*static_cast<double>(n));
  myP = arma::max(mask, myP);
  // myP /= (2.0*static_cast<double>(n));
  // arma::uvec myPsmall = arma::find(myP < 1e-10);
  // myP(myPsmall) = 1e-10;


  // P = (P+P.t());
  // P /= (2*n);
  // P = arma::max(mask,P);

  // 3-2. Main Iteration
  arma::mat Q(n,n,arma::fill::zeros);
  arma::mat PQ(n,n,arma::fill::zeros);
  for (int it=0;it<maxiter;it++){
    dC.zeros();
    Q = computeQ(Y);
    Q = (Q+Q.t())/2;
    Q = arma::max(mask,Q);
    for (int i=0;i<n;i++){
      for (int j=0;j<n;j++){
        if (j!=i){
          dC.col(i) += 4*(myP(i,j)-Q(i,j))*(Y.col(i)-Y.col(j))/(1+pow(norm(Y.col(i)-Y.col(j),2),2));
        }
      }
    }
    y_incs = momentum*y_incs - eta*dC;
    Y += y_incs;
    mat JitMat = jitter*randn<mat>(ndim,n);
    Y += JitMat;
    vec MeanY = mean(Y,1);
    for (int i=0;i<n;i++){
      Y.col(i) -= MeanY;
    }
    jitter *= decay;
  }
  return(Y);
}

// 4. eigenmaps : given weight matrix, compute various embeddings
// [[Rcpp::export]]
Rcpp::List method_eigenmaps(arma::mat& W){
  // 4-1. setting
  const int n = W.n_cols;
  const int m = W.n_rows;
  if (m!=n){
    Rcpp::stop("ERROR : not a square matrix W.");
  }

  // 4-2. compute a normalized graph laplacian
  //    it says, eig_pair is not supported. don't know why.
  vec onesN = ones<vec>(n);
  vec d = W*onesN;
  mat I(n,n,fill::eye);
  mat nL = I - ((diagmat(1/d))*W);

  // 4-3. compute eigenvalues and eigenvectors
  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, nL);
  return Rcpp::List::create(Rcpp::Named("eigval")=eigval,
                            Rcpp::Named("eigvec")=eigvec);
}

// 5. sammon : sammon mapping updates
// X     : (p-by-n) for armadillo convenience
// Yinit : (n-by-d)


// 6. lleW : compute weight matrix W
// [[Rcpp::export]]
arma::vec method_lleW(arma::mat& mat_tgt, arma::vec& vec_tgt, const double regparam){
  // 6-1. basic settings
  const int p = mat_tgt.n_rows;
  const int k = mat_tgt.n_cols;

  // 6-2. compute C
  mat C(k, k, fill::zeros);
  for (int i=0;i<k;i++){
    vec tgtr1 = mat_tgt.col(i);
    vec diff1 = vec_tgt-tgtr1;
    for (int j=i;j<k;j++){
      if (i==j){
        C(i,i) =  dot(diff1,diff1);
      } else {
        vec tgtr2 = mat_tgt.col(j);
        vec diff2 = vec_tgt-tgtr2;
        double elemC = dot(diff1,diff2);
        C(i,j) = elemC;
        C(j,i) = elemC;
      }
    }
  }

  // 6-3. solve for the equation
  vec onesK = ones<vec>(k);
  vec w(k);
  mat I = eye<mat>(k,k);

  if (k>p){
    // 6-3-1. with regularization : k>p
    mat LHS = C + regparam*I;
    w = solve(LHS,onesK);
    w /= sum(w);
  } else {
    // 6-3-2. without regularization
    w = solve(C,onesK);
    w /= sum(w);
  }

  // 6-4. return results
  return(w);
}

// 7. lleWauto : compute weight matrix W with automatic regularization
// [[Rcpp::export]]
Rcpp::List method_lleWauto(arma::mat& mat_tgt, arma::vec& vec_tgt){
  // 7-1. basic settings
  const int p = mat_tgt.n_rows;
  const int k = mat_tgt.n_cols;

  // 7-2. compute C
  mat C(k, k, fill::zeros);
  for (int i=0;i<k;i++){
    vec tgtr1 = mat_tgt.col(i);
    vec diff1 = vec_tgt-tgtr1;
    for (int j=i;j<k;j++){
      if (i==j){
        C(i,i) =  dot(diff1,diff1);
      } else {
        vec tgtr2 = mat_tgt.col(j);
        vec diff2 = vec_tgt-tgtr2;
        double elemC = dot(diff1,diff2);
        C(i,j) = elemC;
        C(j,i) = elemC;
      }
    }
  }

  // 7-3. solve for the equation
  vec onesK = ones<vec>(k);
  vec w(k), wtemp(k);
  mat I(k,k,fill::eye);
  double regparam = 0;
  mat C2 = C.t()*C;

  if (k>p){
    // 7-3-1. with regularization : k>p
    vec tgtlogs   = logspace<vec>(-2,2,10);
    double gcvval = 123456789;
    for (int it=0;it<10;it++){
      double regtmp = tgtlogs(it);
      mat LHS = C + regtmp*I;
      wtemp = solve(LHS,onesK);

      // main GCF computation
      double nominator   = pow(norm(C*wtemp-onesK),2);
      double denominator = pow(trace(I-C*(solve(C2+regtmp*I,C.t()))),2);
      double gcvtmp = nominator/denominator;
      if (gcvtmp <= gcvval){
        gcvval   = gcvtmp;
        regparam = regtmp;
        w        = wtemp;
        w       /= sum(w);
      }
    }
  } else {
    // 7-3-2. without regularization
    w = solve(C,onesK);
    w /= sum(w);
  }

  // 7-4. return results
  return Rcpp::List::create(Rcpp::Named("w")=w,
                            Rcpp::Named("regparam")=regparam);
}


// 8. lleM : main 2 for computing low-D embedding
// [[Rcpp::export]]
Rcpp::List method_lleM(arma::mat& W){
  const int n = W.n_cols;
  mat I(n,n,fill::eye);
  mat M = I.t()*I - W.t()*I - I*W + W.t()*W;

  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, M);
  return Rcpp::List::create(Rcpp::Named("eigval")=eigval,
                            Rcpp::Named("eigvec")=eigvec);
}


// 9. REE : Robust Euclidean Embedding
arma::mat method_ree_subgradient(arma::mat B, arma::mat W, arma::mat D){
  const int n = B.n_cols;
  arma::mat GB(n,n,fill::zeros);
  arma::mat distB(n,n,fill::zeros);
  // fill in the distB
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      distB(i,j) = B(i,i)+B(j,j)-B(i,j)-B(j,i);
    }
  }
  // main computation
  // off diagonals first
  double tgt = 0.0;
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      if (i!=j){
        tgt = distB(i,j);
        if (tgt<D(i,j)){
          GB(i,j) = W(i,j)*1.0;
        } else {
          GB(i,j) = W(i,j)*(-1.0);
        }
      }
    }
  }
  // diagonals
  for (int i=0;i<n;i++){
    tgt = 0.0;
    for (int k=0;k<n;k++){
      if (distB(i,k)>D(i,k)){
        tgt += W(i,k);
      } else {
        tgt -= W(i,k);
      }
    }
    GB(i,i) = tgt;
  }
  return(GB);
}

double method_ree_cost(arma::mat W, arma::mat D, arma::mat B){
  const int n = W.n_cols;
  double score = 0.0;
  double term1 = 0.0;
  double term2 = 0.0;
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      term1 = W(i,j);
      term2 = std::abs(D(i,j)-(B(i,i)+B(j,j)-B(i,j)-B(j,i)));
      score += (term1*term2);
    }
  }
  return(score);
}

// [[Rcpp::export]]
Rcpp::List method_ree(arma::mat& B, arma::mat& W, arma::mat& D, const double initc,
                      const double abstol, const int maxiter){
  // 1. settings
  const int n = B.n_cols;
  arma::mat Bold = B;
  arma::mat Bnew(n,n,fill::zeros);
  arma::mat Btmp(n,n,fill::zeros);
  arma::mat Btgt(n,n,fill::zeros);

  // 2. iterate
  double alpha  = 0.0;  // subgradient stepsize
  double gap    = 10.0;
  int    iter   = 1;
  double cost   = 10000.0;
  double cost_k = 0.0;
  while (gap > abstol){
    // 2-1. update
    alpha = initc/static_cast<double>(iter);
    Btmp  = Bold - alpha*method_ree_subgradient(Bold, W, D);

    // 2-2. spectral decomposition
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, Btmp);

    // 2-3. eigenvalue removal of zeros
    for (int i=0;i<n;i++){
      if (eigval(i)<0){
        eigval(i)=0;
      }
    }
    Bnew = eigvec*diagmat(eigval)*eigvec.t();

    // 2-4. compute newer cost
    cost_k = method_ree_cost(W,D,Bnew);
    if (cost_k < cost){
      cost = cost_k;
      Btgt = Bnew;
    }

    // 2-5. stopping criterion and iteration update
    iter += 1;
    gap  = norm(Bnew-Bold,"f");
    Bold = Bnew;
    if (iter >= maxiter){
      gap = abstol/10.0;
    }
  }

  // 3. return results
  return Rcpp::List::create(Rcpp::Named("B")=Bold,
                            Rcpp::Named("iter")=iter);
}


// 10. SPE : Stochastic Proximity Embedding
// [[Rcpp::export]]
arma::mat method_spe(arma::mat& R, arma::mat& iX, const int C, const int S,
                     double lambda, double drate, arma::mat matselector){
  // 1. setup
  arma::mat X = iX;
  // const int nrowX = iX.n_rows; // unused flag
  const double tolerance = 0.00001;
  double denomeps = 0.0000001;

  // 2. let's iterate
  arma::vec idx2;
  int selectit = 0;
  int i=0;
  int j=0;
  double dij = 0.0;
  arma::rowvec xi;
  arma::rowvec xj;
  double coefficient = 0.0;
  for (int iterC=0;iterC<C;iterC++){ // iteration for 'C' cycles
    // 2-1. iteration for S steps
    for (int iterS=0;iterS<S;iterS++){
      // 2-1-1. select two points
      i = static_cast<int>(matselector(selectit, 0));
      j = static_cast<int>(matselector(selectit, 1));

      xi = X.row(i);
      xj = X.row(j);

      // 2-1-2. compute distance dij
      dij = arma::norm(xi-xj, 2);

      // 2-1-3. maybe, update?
      if (std::abs(dij-R(i,j)) > tolerance){
        coefficient = (lambda*(R(i,j)-dij))/(2*(dij+denomeps));

        X.row(i) = xi + coefficient*(xi-xj);
        X.row(j) = xj + coefficient*(xj-xi);
      }
      // 2-1-10. update selectit
      selectit += 1;
    }
    // 2-2. update lambda
    lambda *= drate;
  }

  // 3. return output
  return(X);
}
// [[Rcpp::export]]
arma::mat method_ispe(arma::mat& R, arma::mat& iX, const int C, const int S,
                      double lambda, double drate, arma::mat matselector, const double cutoff){
  // 1. setup
  arma::mat X = iX;
  // const int nrowX = iX.n_rows;      // unused flag on windows
  // const double tolerance = 0.00001; // unused flag
  double denomeps = 0.0000001;

  // 2. let's iterate
  arma::vec idx2;
  int selectit = 0;
  int i=0;
  int j=0;
  double dij = 0.0;
  arma::rowvec xi;
  arma::rowvec xj;
  double coefficient = 0.0;
  for (int iterC=0;iterC<C;iterC++){ // iteration for 'C' cycles
    // 2-1. iteration for S steps
    for (int iterS=0;iterS<S;iterS++){
      // 2-1-1. select two points
      i = static_cast<int>(matselector(selectit, 0));
      j = static_cast<int>(matselector(selectit, 1));

      xi = X.row(i);
      xj = X.row(j);

      // 2-1-2. compute distance dij
      dij = arma::norm(xi-xj, 2);

      // 2-1-3. maybe, update?
      if ((R(i,j)<=cutoff)||(dij<R(i,j))){
        coefficient = (lambda*(R(i,j)-dij))/(2*(dij+denomeps));
        X.row(i) = xi + coefficient*(xi-xj);
        X.row(j) = xj + coefficient*(xj-xi);
      }

      // 2-1-10. update selectit
      selectit += 1;
    }
    // 2-2. update lambda
    lambda *= drate;
  }

  // 3. return output
  return(X);
}

// 11. CRCA : Curvilinear Component Analysis
arma::mat method_crca_dist(arma::mat RowMat){
  const int n = RowMat.n_rows;
  arma::mat output(n,n,fill::zeros);
  arma::rowvec vec1;
  arma::rowvec vec2;
  double norm12 = 0.0;
  for (int i=0;i<(n-1);i++){
    vec1 = RowMat.row(i);
    for (int j=(i+1);j<n;j++){
      vec2 = RowMat.row(j);

      norm12 = arma::norm(vec1-vec2,2);
      output(i,j) = norm12;
      output(j,i) = norm12;
    }
  }
  return(output);
}
// [[Rcpp::export]]
Rcpp::List method_crca(arma::mat& Xij, arma::mat& Yinit, double lambda, double alpha, const int maxiter, const double tolerance, arma::vec& vecselector){
  // 1. get parameters
  const int n = Yinit.n_rows;
  // const int ndim = Yinit.n_cols; :: unused flag from windows

  // 2. settings
  arma::mat Y = Yinit; // deep copy
  arma::mat Yij = method_crca_dist(Y);

  // 3. iterate !
  double increment = 1000.0;
  int t = 0;
  double tdb = 0.0;
  int i = 0;
  double alpha_t = 0.0;

  arma::rowvec veci;
  arma::rowvec vecj;

  while (increment > tolerance){
    // 3-1. select i
    i    = static_cast<int>(vecselector(t));
    veci = Y.row(i);

    // 3-2. alpha_t
    tdb = static_cast<double>(t);
    alpha_t = alpha/(1.0+tdb);

    // 3-3. iterate over j
    for (int j=0;j<n;j++){
      vecj = Y.row(j);
      if (i!=j){
        if (Yij(i,j)<=lambda){
          Y.row(j) = vecj + alpha_t*(Xij(i,j)-Yij(i,j))*(vecj-veci)/(Yij(i,j));
        }
      }
    }

    // 3-4. update Yij and increment //////////////////////////////////////////// abs part
    Yij = method_crca_dist(Y);
    increment = (abs(Xij-Yij)).max();

    // 3-5. update iteration count
    t = t + 1;
    if (t>=maxiter){
      break;
    }
  }

  // 4. return output
  return Rcpp::List::create(Rcpp::Named("Y")=Y,
                            Rcpp::Named("niter")=t);
}

// 12. BMDS : Bayesian MDS
// [[Rcpp::export]]
double bmds_compute_SSR(arma::mat &D, arma::mat &Delta){
  // parameters
  int N = D.n_rows;
  double NN = static_cast<double>(N);

  // compute via iteration
  double outval = 0.0;
  double tobesq = 0.0;
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      tobesq = (D(i,j)-Delta(i,j));
      outval += (tobesq*tobesq)/NN;
    }
  }
  return(outval);
}
double bmds_compute_SSR_xmat(arma::mat &D, arma::mat &Xnew){ // this one is using matrix data
  int N = D.n_rows; double NN = static_cast<double>(N);
  int p = Xnew.n_cols;

  double outval = 0.0;
  double tobesq = 0.0;

  arma::rowvec xvec1(p,fill::zeros);
  arma::rowvec xvec2(p,fill::zeros);

  double Delij = 0.0;
  for (int i=0;i<N;i++){
    xvec1 = Xnew.row(i);
    for (int j=(i+1);j<N;j++){
      xvec2  = Xnew.row(j);
      Delij  = arma::norm(xvec1-xvec2, 2);
      tobesq = D(i,j)-Delij;
      outval+= (tobesq*tobesq)/NN;
    }
  }
  return(outval);
}
arma::mat bmds_compute_pdmat(arma::mat &X){
  int N = X.n_rows;
  int p = X.n_cols;
  arma::mat output(N,N,fill::zeros);
  arma::vec tgt1(p,fill::zeros);
  arma::vec tgt2(p,fill::zeros);
  double tmpval = 0.0;
  for (int i=0;i<(N-1);i++){
    tgt1 = X.row(i).t();
    for (int j=0;j<N;j++){
      tgt2 = X.row(j).t();
      tmpval = arma::norm(tgt1-tgt2,2);
      output(i,j) = tmpval;
      output(j,i) = tmpval;
    }
  }
  return(output);
}
arma::mat bmds_crotX(arma::mat X){
  int N = X.n_rows;
  int p = X.n_cols;

  arma::mat Xtmp(N,p,fill::zeros);
  arma::rowvec xmean = arma::mean(X, 0);
  for (int i=0;i<N;i++){
    Xtmp.row(i) = X.row(i)-xmean;
  }

  arma::mat Xcov = Xtmp.t()*Xtmp/(static_cast<double>(N));

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Xcov);

  arma::mat output = Xtmp*eigvec;
  return(output);
}
arma::rowvec bmds_update_xvec(arma::mat D, arma::mat X, int id, double sigma2, double constant, arma::mat Lbdmat){
  int N = X.n_rows; double NN = static_cast<double>(N);
  int p = X.n_cols;
  arma::mat Xold = X;
  arma::mat Xtgt = X;
  double stepsize = static_cast<double>(std::sqrt(static_cast<float>(sigma2*constant/(NN-1.0))));
  for (int i=0;i<p;i++){
    Xtgt(id,i) += R::rnorm(0.0, stepsize);
  }
  double sigma = std::sqrt(static_cast<float>(sigma2));

  arma::vec xtgt = Xtgt.row(id).t(); // column vectors
  arma::vec xold = Xold.row(id).t();

  // common variables
  double tmpval = 0.0;

  // need to evaluate two ratio
  // (1) compute for xtgt
  arma::mat Deltgt = bmds_compute_pdmat(Xtgt);
  double Q1tgt = 0.0;
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      tmpval = D(i,j)-Deltgt(i,j);
      Q1tgt += (tmpval*tmpval)/sigma2;
    }
  }
  double Q2tgt = arma::dot(xtgt, arma::solve(Lbdmat, xtgt));
  double t3tgt = 0.0;
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      if (i!=j){
        t3tgt += static_cast<double>(std::sqrt(static_cast<float>(R::pnorm5(Deltgt(i,j)/sigma,0.0,1.0,1,0))));
      }
    }
  }
  double ftgt = -(Q1tgt+Q2tgt)/2.0 - t3tgt;

  // (2) compute for xold
  arma::mat Delold = bmds_compute_pdmat(Xold);
  double Q1old = 0.0;
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      tmpval = D(i,j)-Delold(i,j);
      Q1old += (tmpval*tmpval)/sigma2;
    }
  }
  double Q2old = arma::dot(xold, arma::solve(Lbdmat, xold));
  double t3old = 0.0;
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      if (i!=j){
        t3old += static_cast<double>(std::sqrt(static_cast<float>(R::pnorm5(Delold(i,j)/sigma,0.0,1.0,1,0))));
      }
    }
  }
  double fold = -(Q1old+Q2old)/2.0 - t3old;

  // (3) compute the ratio (?)
  double fratio = std::exp(static_cast<float>(ftgt-fold));
  if (fratio >= 1){
    fratio = 1.0;
  }
  double rnumbr = R::runif(0.0, 1.0);
  if (rnumbr <= fratio){ // accept
    return(xtgt.t());
  } else {
    return(xold.t());
  }
}
double my_invgamma(double alpha, double beta){
  return(1.0/R::rgamma(alpha,1.0/beta));
}
double my_dinvgamma(double x, double alpha, double beta){
  return(1.0/R::dgamma(x, alpha, 1.0/beta, 0));
}
// [[Rcpp::export]]
Rcpp::List main_bmds(arma::mat D, arma::mat X0, double sigg0,
                     double a, double alpha, int maxiter, double constant, bool verbose,
                     arma::vec betas){
  // 1) some parameters
  int N = X0.n_rows; double NN = static_cast<double>(N);
  int p = X0.n_cols;
  double m = NN*(NN-1.0)/2.0;


  // 2) setup
  arma::mat Xold = bmds_crotX(X0); // X will not be recorded, just use
  arma::mat Xnew(N,p,fill::zeros);
  arma::mat Xsol = Xold;

  double SSRnew = 0.0;
  double SSRold = bmds_compute_SSR_xmat(D, Xold);
  double SSRsol = SSRold;

  arma::mat Sold(p,p,fill::zeros);
  double sigma2 = sigg0;
  double sigtmp = 0.0;
  arma::vec vecs(p,fill::zeros);
  arma::vec lambdas(p,fill::zeros);
  arma::mat Lbdmat;
  arma::rowvec tmprow(p,fill::zeros);
  double b = (a-1)*SSRold/m; // paper's setup

  double varalpha = 0.0;
  double varbeta  = 0.0;
  double varvar = 0.0;
  double varratio = 0.0;

  // 3) iteration
  int accept = 0;
  for (int i=0;i<maxiter;i++){
    // 3-1. update lambdas
    for (int j=0;j<p;j++){ // compute sample variances for each coordinate
      vecs(j) = arma::var(Xold.col(j))*NN;
    }
    Sold = Xold.t()*Xold/NN;
    for (int j=0;j<p;j++){ // sample from IG
      lambdas(j) = my_invgamma(alpha+NN/2.0, betas(j) + vecs(j)/2.0); // according to the paper's choice
    }
    Lbdmat = arma::diagmat(lambdas);

    // 3-2. update X
    Xnew   = Xold;
    for (int j=0;j<N;j++){ // for each row
      tmprow = bmds_update_xvec(D, Xnew, j, sigma2, constant, Lbdmat);
      Xnew.row(j) = tmprow;
    }
    SSRnew = bmds_compute_SSR_xmat(D, Xnew); // update SSR
    Xnew   = bmds_crotX(Xnew); // centering + rotation

    // 3-3. update sigma using MH
    varalpha = m/2 + a;
    varbeta  = SSRnew/2 + b;
    varvar   = (varbeta*varbeta)/((varalpha-1)*(varalpha-1)*(varalpha-2));

    sigtmp = sigma2 + R::rnorm(0, static_cast<double>(std::sqrt(static_cast<float>(constant*varvar))));
    if (sigtmp > 0){ // let's compare
      varratio = my_dinvgamma(sigtmp,varalpha,varbeta)/my_dinvgamma(sigma2,varalpha,varbeta);
      if (varratio > 1){
        varratio = 1.0;
      }
      if (R::runif(0,1) <= varratio){
        sigma2 = sigtmp;
      }
    }


    // 3-4. update correspondingly
    if (SSRnew < SSRsol){ // running record of the best solution
      SSRsol = SSRnew;
      Xsol   = Xnew;
    }
    SSRold = SSRnew;
    Xold   = Xnew;

    // 3-5. report the update
    if (verbose==true){
      Rcpp::Rcout << "** do.bmds : iteration " << i+1 << "/" << maxiter << " complete." << std::endl;
    }
  }

  // 4) return
  return Rcpp::List::create(Rcpp::Named("solX")=Xsol);
}
