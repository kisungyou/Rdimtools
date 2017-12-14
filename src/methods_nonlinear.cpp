#include <RcppArmadillo.h>
#include "methods_nonlinear.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// subroutines
//    1. Compute Q for Low-Dimensional Values
arma::mat computeQ(arma::mat& Y){
  const int n = Y.n_cols;
  mat Q(n,n);
  for (int i=0;i<n;i++){
    double cdenom = 0;
    for (int k=0;k<n;k++){
      if (i==k){
        cdenom += 0;
      } else{
        cdenom += exp(-pow(norm(Y.col(i)-Y.col(k)),2));
      }
    }
    for (int j=0;j<n;j++){
      if (i==j){
        Q(i,j) = 0;
      } else {
        Q(i,j) = exp(-pow(norm(Y.col(i)-Y.col(j)),2))/cdenom;
      }
    }
  }
  return(Q);
}

// 1. SNE : Stochastic Neighbor Embedding
//' @keywords internal
// [[Rcpp::export]]
arma::mat method_sne(arma::mat& P, const int ndim, const double eta,
                     const int maxiter, double jitter, double decay,
                     const double momentum){
  // 1-1. Initialize
  const int n = P.n_cols;
  mat Y = 0.0001*randn<mat>(ndim,n);
  mat dC(ndim,n);
  mat y_incs(ndim,n);
  mat mask = (1e-10)*ones<mat>(n,n);
  P = arma::max(mask,P);


  // 1-2. Main Iteration
  mat Q(n,n);
  mat PQ(n,n);
  for (int it=0;it<maxiter;it++){
    dC.zeros();
    Q = computeQ(Y);
    Q = arma::max(mask,Q);
    for (int i=0;i<n;i++){
      for (int j=0;j<n;j++){
        if (j!=i){
          dC.col(i) += 2*(P(j,i)+P(i,j)-Q(i,j)-Q(j,i))*(Y.col(i)-Y.col(j));
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


// 2. Symmetric SNE : Stochastic Neighbor Embedding
//' @keywords internal
// [[Rcpp::export]]
arma::mat method_snesym(arma::mat& P, const int ndim, const double eta,
                     const int maxiter, double jitter, double decay,
                     const double momentum){
  // (1) Initialize
  const int n = P.n_cols;
  mat Y = 0.0001*randn<mat>(ndim,n);
  mat dC(ndim,n);
  mat y_incs(ndim,n);
  mat mask = (1e-10)*ones<mat>(n,n);
  P = (P+P.t());
  P /= (2*n);
  P = arma::max(mask,P);

  // (2) Main Iteration
  mat Q(n,n);
  mat PQ(n,n);
  for (int it=0;it<maxiter;it++){
    dC.zeros();
    Q = computeQ(Y);
    Q = (Q+Q.t())/2;
    Q = arma::max(mask,Q);
    for (int i=0;i<n;i++){
      for (int j=0;j<n;j++){
        if (j!=i){
          dC.col(i) += 4*(P(i,j)-Q(i,j))*(Y.col(i)-Y.col(j));
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

// 3. tSNE : t-Stochastic Neighbor Embedding
//' @keywords internal
// [[Rcpp::export]]
arma::mat method_tsne(arma::mat& P, const int ndim, const double eta,
                        const int maxiter, double jitter, double decay,
                        const double momentum){
  // 3-1. Initialize
  const int n = P.n_cols;
  mat Y = 0.0001*randn<mat>(ndim,n);
  mat dC(ndim,n);
  mat y_incs(ndim,n);
  mat mask = (1e-10)*ones<mat>(n,n);
  P = (P+P.t());
  P /= (2*n);
  P = arma::max(mask,P);

  // 3-2. Main Iteration
  mat Q(n,n);
  mat PQ(n,n);
  for (int it=0;it<maxiter;it++){
    dC.zeros();
    Q = computeQ(Y);
    Q = (Q+Q.t())/2;
    Q = arma::max(mask,Q);
    for (int i=0;i<n;i++){
      for (int j=0;j<n;j++){
        if (j!=i){
          dC.col(i) += 4*(P(i,j)-Q(i,j))*(Y.col(i)-Y.col(j))/(1+pow(norm(Y.col(i)-Y.col(j),2),2));
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
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_eigenmaps(arma::mat& W){
  // 4-1. setting
  const int n = W.n_cols;
  if (W.n_rows!=n){
    Rcpp::stop("ERROR : not a symmetric one in size.");
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
//' @keywords internal
// [[Rcpp::export]]
arma::mat method_sammon(arma::mat& X, arma::mat& Yinit){
  // 5-1. basic settings
  const int n = Yinit.n_rows;
  const int d = Yinit.n_cols;

  // 5-2. iteration prep
  mat Yold = Yinit;
  mat Ynew(n,d);

  const int maxiter = 1000;
  const double thrstop = 1e-10;
  double incnorm;

  // 5-3. common values from original matrix X
  mat DPJX(n,n);
  double c = 0;
  for (int i=0;i<n;i++){
    for (int j=i;j<n;j++){
      if (i!=j){
        double DPJXstar = norm(X.col(i)-X.col(j),2);
        DPJX(i,j) = DPJXstar;
        DPJX(j,i) = DPJXstar;
        c += DPJXstar;
      }
    }
  }

  // 5-4. main iteration
  mat DPJY(n,n);
  for (int m=0;m<maxiter;m++){
   //  5-4-1. define & compute new DPJ for current iterate
    DPJY.zeros();
    for (int i=0;i<n;i++){
      for (int j=i;j<n;j++){
        if (i!=j){
          double DPJcurrent = norm(Yold.row(i)-Yold.row(j));
          DPJY(i,j) = DPJcurrent;
          DPJY(j,i) = DPJcurrent;
        }
      }
    }
    //  5-4-2. update with for Ynew
    for (int p=0;p<n;p++){
      for (int q=0;q<d;q++){
        double der1 = 0;
        double der2 = 0;

        for (int j=0;j<n;j++){
          if (p!=j){
            double dpjstar = DPJX(p,j);
            double dpj     = DPJY(p,j);
            der1 += ((dpjstar-dpj)/(dpjstar*dpj))*(Yold(p,q)-Yold(j,q));
            der2 += (dpjstar-dpj)-(pow((Yold(p,q)-Yold(j,q)),2)/dpj)*(1+((dpjstar-dpj)/dpj));
          }
        }

        der1 *= (-2/c);
        der2 *= (-2/c);
        if (der2 < 0){
          der2 *= -1;
        }

        Ynew(p,q) = Yold(p,q) - 0.3*(der1/der2);
      }
    }

    //  5-4-3. Break with incremental norm condition
    incnorm = norm(Yold-Ynew);
    if (incnorm < thrstop){
      Yold = Ynew;
      break;
    }

    //  5-4-4. update
    Yold = Ynew;
  }

  // 5-5. return output
  return(Yold);
}



// 6. lleW : compute weight matrix W
//' @keywords internal
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
//' @keywords internal
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
//' @keywords internal
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
//' @keywords internal
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

//' @keywords internal
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

//' @keywords internal
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
      cost = cost;
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
  const int nrowX = iX.n_rows;
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
  const int nrowX = iX.n_rows;
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
