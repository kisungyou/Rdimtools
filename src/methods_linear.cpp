/*
 * 01. PCA
 * 02. MDS
 * 03. MDS given D
 * 04. ICA
 * 05. RNDPROJ
 * 06. FA
 * 07. LPP
 * 08. NPE
 * 09. OLPP
 * 10. BPCA
 * 11. EXTLPP
 * 12. LSPP
 * 13. KMCC
 * 14. LFDA
 * 15. NNPROJMAX & NNPROJMIN
 * 16. NNEMBEDMIN
 * 17. SPUFS
*/

#include <RcppArmadillo.h>
#include "methods_linear.h"
#include "methods_handytools.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// 1. PCA : Principal Component Analysis
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_pca(arma::mat& psdX){
  // define auxiliary stuffs
  arma::vec eigval;
  arma::mat eigvec;
  // run eigendecomposition
  eig_sym(eigval, eigvec, psdX);
  // return output
  return Rcpp::List::create(Rcpp::Named("eigval")=eigval,
                            Rcpp::Named("eigvec")=eigvec);
}

// 2. MDS : Multi-dimensional Scaling (Classical)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_mds(arma::mat& centerX){
  // const int n = centerX.n_cols; // unused flag
  // X(centered) = U*S*V.t()
  mat U;
  vec s;
  mat V;
  svd(U, s, V, centerX);
  // output
  return Rcpp::List::create(Rcpp::Named("eigval")=s,
                            Rcpp::Named("eigvec")=V);
}

// 3. MDS given D
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_mdsD(arma::mat& D){
  const int n = D.n_cols;
  // 3-1. H = I-(1/n)*ones(n,n)
  mat D2 = pow(D,2);
  mat I(n,n); I.eye();
  mat oneN(n,n); oneN.ones(); oneN /= n;
  mat H = I - oneN;

  // 3-2. B
  mat B = (-0.5)*(H*D2*H.t());

  // 3-3. eigenvalue decomposition
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, B);

  // 3-4. Return !
  return Rcpp::List::create(Rcpp::Named("eigval")=eigval,
                            Rcpp::Named("eigvec")=eigvec);
}


// 4. ICA : (Fast) Independent Component Analysis
// tnum = 1 (logcosh), 2 (exp), 3 (poly)
// sym  = FALSE (just run) / TRUE (decorrelate)

typedef vec (*icaPtr)(const vec& x, const double tpar);
vec ica_logcosh(const vec& x, const double tpar){
  vec y = tanh(tpar*x);
  return(y);
}
vec ica_logcoshp(const vec& x, const double tpar){
  vec y = tpar*(1-pow(tanh(tpar*x),2));
  return(y);
}
vec ica_exp(const vec& x, const double tpar){
  const int n = x.n_elem;
  vec y(n,fill::zeros);
  for (int i=0;i<n;i++){
    y(i) = x(i)*exp(-tpar*pow(x(i),2)/2);
  }
  //vec y =x*exp(-tpar*pow(x,2)/2);
  return(y);
}
vec ica_expp(const vec& x, const double tpar){
  const int n = x.n_elem;
  vec y(n,fill::zeros);
  for (int i=0;i<n;i++){
    double x2 = pow(x(i),2);
    y(i) = (1-tpar*x2)*exp((-tpar*x2)/2);
  }
  /*vec ones;// exp is the problem
  ones.copy_size(x);
  ones.ones();
  vec x2 = pow(x,2);
  vec y = (ones-tpar*x2)*exp((-tpar*x2)/2);*/
  return(y);
}
vec ica_poly(const vec& x, const double tpar){
  vec y = pow(x,3);
  return(y);
}
vec ica_polyp(const vec& x, const double tpar){
  vec y = 3*(pow(x,2));
  return(y);
}

XPtr<icaPtr> decideICAg(const int n){
  if (n==1){
    return(XPtr<icaPtr>(new icaPtr(&ica_logcosh)));
  } else if (n==2){
    return(XPtr<icaPtr>(new icaPtr(&ica_exp)));
  } else if (n==3){
    return(XPtr<icaPtr>(new icaPtr(&ica_poly)));
  } else {
    return XPtr<icaPtr>(R_NilValue);
  }
}
XPtr<icaPtr> decideICAgprime(const int n){
  if (n==1){
    return(XPtr<icaPtr>(new icaPtr(&ica_logcoshp)));
  } else if (n==2){
    return(XPtr<icaPtr>(new icaPtr(&ica_expp)));
  } else if (n==3){
    return(XPtr<icaPtr>(new icaPtr(&ica_polyp)));
  } else {
    return XPtr<icaPtr>(R_NilValue);
  }
}

//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_ica(arma::mat& X, const int C, const int maxiter,
                      const double tol, const int tnum, const double tpar,
                      bool sym){
  // 1. setting
  const int N = X.n_rows;
  const int M = X.n_cols;
  mat W       = 0.01*(randu<mat>(N,C));
  for (int i=0;i<C;i++){
    W.col(i) /= norm(W.col(i),2);
  }

  // 2. use function pointer for case branching
  XPtr<icaPtr> xpfun1 = decideICAg(tnum);
  XPtr<icaPtr> xpfun2 = decideICAgprime(tnum);

  icaPtr trf_g  = *xpfun1; // transformation for g
  icaPtr trf_gp = *xpfun2; // transformation for gprime

  vec onesM(M);
  onesM.ones();

  // 3. main computation
  for (int it=0;it<C;it++){
    vec wold = W.col(it);
    vec wnew, wtX(N), wtX1(N), wtX2(N);
    double tgtgap  = 100.0;
    int    tgtiter = 0;

    while (tgtgap > tol){
      wtX = (X.t()*wold);
      wtX1 = trf_g(wtX,tpar);
      wtX2 = trf_gp(wtX,tpar);

      wnew = (X*(wtX1))/M - as_scalar((wtX2.t()*onesM))*wold/M;
      //wnew = randu<vec>(N);
      if (it>0){
        for (int i=0;i<(it-1);i++){
          wnew -= as_scalar((wnew.t())*W.col(i))*W.col(i);
        }
      }
      wnew /= norm(wnew);

      tgtiter += 1;
      tgtgap = norm(wnew-wold);

      wold = wnew;
      if (tgtiter>maxiter){
        break;
      }
    }
    W.col(it) = wold;
  }

  // 4. symmetric decorrelation
  if (sym==false){
    mat S = (W.t())*X;
    return Rcpp::List::create(Rcpp::Named("S")=S,
                              Rcpp::Named("W")=W);
  } else if (sym==true){
    mat WWT = W*W.t();
    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, WWT);
    mat preW = eigvec*diagmat(1/sqrt(eigval))*eigvec.t();
    mat Wnew = preW*W;

    mat S = (Wnew.t())*X;
    return Rcpp::List::create(Rcpp::Named("S")=S,
                              Rcpp::Named("W")=Wnew);
  } else {
    Rcpp::stop(" fastica : sym flag is incorrect.");
  }
}


// 5. rpgauss : Random Projection for Gaussian Case
//    Note that in this case, I simply followed R manner of data arrangement
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_rpgauss(arma::mat& X, const int k){
  // 5-1. setting
  // const int n = X.n_rows; // unused flag
  const int d = X.n_cols;
  mat R(d,k);
  mat G = (randu<mat>(d,k));

  // 5-2. generate and iterate
  R.col(0) = G.col(0)/norm(G.col(0));

  for (int i=1;i<k;i++){
    vec tgt = G.col(i);
    for (int j=0;j<i;j++){
      vec u = R.col(j);
      tgt -= (dot(u,tgt)/dot(tgt,tgt))*u;
    }
    R.col(i) = tgt/norm(tgt);
  }

  // 5-3. Compute
  mat Y = X*R;
  return Rcpp::List::create(Rcpp::Named("R")=R,
                            Rcpp::Named("Y")=Y);
}

/*
 * 6. Factor Analysis
 *    Notations are consistent with that of Gharahmani's paper.
 */
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_fa(arma::mat& X, const int k, const int maxiter, const double tolerance){
  // 6-1. basic settings
  const int p = X.n_rows;
  const int n = X.n_cols;

  // 6-2. other defined materials
  mat Ez(k,n);
  mat beta(k,p);
  mat LPhi(p,p);
  mat eyeK(k,k,fill::eye);
  mat Pinv(p,p,fill::zeros);
  mat EzzSum(k,k);
  double inctol = 0;

  // 6-3. initialization
  mat Lold = randn<mat>(p,k);
  mat Lnew(p,k,fill::zeros);
  vec Pold = ones<vec>(p);
  mat Ptmp(p,p,fill::zeros);
  vec Pnew(p);

  // 6-4. main iteration
  for (int it=0;it<maxiter;it++){
    // 6-4-E1. common again
    Pinv = diagmat(1/Pold);
    // 6-4-E2. inverse of LL^T+\Phi : LPhi
    LPhi = Pinv - Pinv*Lold*solve(eyeK+Lold.t()*Pinv*Lold,Lold.t()*Pinv);
    // 6-4-E3. beta
    beta = Lold.t()*LPhi;
    // 6-4-E4. EZ = [EZ(x1), EZ(x2), ... , EZ(xn)]
    Ez = beta*X;
    // 6-4-E5. EZZsum = sum(EZZ(xi))
    EzzSum = n*(eyeK-beta*Lold) + beta*X*X.t()*beta.t();

    // 6-4-M1. update Lambda : Lnew
    Lnew = (X*Ez.t())*pinv(EzzSum);
    // 6-4-M2. update Phi    : Pnew from Ptmp
    Ptmp = (X - Lnew*Ez)*X.t()/n;
    Pnew = Ptmp.diag();

    // 6-4-Update
    inctol = norm(Lnew-Lold,"fro");
    Lold = Lnew;
    Pold = Pnew;
    if (inctol < tolerance){
      break;
    }
  }

  // 6-5. MLE solution for Z
  mat Z(k,n);
  Pinv = diagmat(1/Pold);
  Z = solve(Lold.t()*Pinv*Lold,Lold.t()*sqrt(Pinv)*X);

  // 6-6. return results
  return Rcpp::List::create(Rcpp::Named("L")=Lold,
                            Rcpp::Named("Z")=Z,
                            Rcpp::Named("Pvec")=Pold);
}

/*
 * 07. Locality Preserving Projections (LPP)
 */

/*
 * 08. Neighborhood Preserving Embedding (NPE)
 */
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_npe(arma::mat& X, arma::mat& W){
  // 08-1. basic settings
  // const int n = X.n_rows; // unused flag
  const int m = X.n_cols;
  const int w1 = W.n_rows;
  const int w2 = W.n_cols;

  if (w1!=w2){
    Rcpp::stop("ERROR : W is not a square matrix.");
  }
  if (m!=w2){
    Rcpp::stop("ERROR : two inputs are not matching.");
  }

  // 08-2. preliminary items
  mat Im(m,m,fill::eye);
  mat M = ((Im-W).t())*(Im-W);

  // 08-3. LHS, RHS and SOL
  mat LHS = X*M*X.t();
  mat RHS = X*X.t();
  mat SOL = solve(RHS, LHS);

  // 08-4. main computation
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, SOL);

  // 08-5. return results
  return Rcpp::List::create(Rcpp::Named("eigval")=eigval,
                            Rcpp::Named("eigvec")=eigvec);
}

/*
 * 09. Orthogonal Locality Preserving Projection (OLPP)
 * NOTE : X in (p-by-n) type for consistency to the notes.
 */
//' @keywords internal
// [[Rcpp::export]]
arma::mat method_olpp(arma::mat& X, arma::mat& S, const int ndim){
  // 9-1. basic settings
  // const int n = S.n_cols; // unused flag
  const int p = X.n_rows;
  // 9-2. laplacian
  arma::mat D = arma::diagmat(sum(S,1));
  arma::mat L = D-S;
  // 9-3. complementary setups
  arma::mat XDXT = X*D*X.t();
  arma::mat XLXT = X*L*X.t();
  // 9-4. ready for accumulations
  arma::mat A(p,ndim,fill::zeros);
  // 9-5. compute the first eigenvector
  arma::vec vals1;
  arma::mat vecs1;
  arma::mat XDLXT = arma::solve(XDXT, XLXT);
  arma::eig_sym(vals1, vecs1, XDLXT);
  A.col(0) = vecs1.col(0);
  // 9-6. main iterations
  arma::mat Ak1; // this is for A(k-1)
  arma::mat Bk1;
  arma::mat Mk;
  arma::mat eyeP(p,p,fill::eye);
  arma::vec valsi; // eigendecomposition of Mk
  arma::mat vecsi; // eigendecomposition of Mk
  for (int i=1;i<ndim;i++){
    // 9-6-1. subsetting
    Ak1 = A.cols(0,i);
    Bk1 = Ak1.t()*(arma::solve(XDXT, Ak1));
    // 9-6-2. compute Mk
    Mk = (eyeP - arma::solve(XDXT, Ak1)*arma::solve(Bk1, Ak1.t()))*XDLXT;
    // 9-6-3. get the smallest eigenvector
    eig_sym(valsi, vecsi, Mk);
    A.col(i) = vecsi.col(0);
  }

  // 9-7. return results
  return(A);
}

/*
 * 10. Bayesian Principal Component Analysis (BPCA)
 * NOTE : T in (d-by-N) type for consistency to the notes.
 */
arma::mat auxiliary_outer(arma::colvec x, arma::colvec y){
  const int nx = x.n_elem;
  const int ny = y.n_elem;
  arma::mat output(nx,ny,fill::zeros);
  for (int i=0;i<ny;i++){
    output.col(i) = x*arma::as_scalar(y(i));
  }
  return(output);
}
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_bpca(arma::mat& T, const double reltol, const int maxiter){
  // 10-1. settings and preliminaries
  int N = T.n_cols;
  int d = T.n_rows;
  int q = d-1;

  // 10-2. initialize
  arma::mat Wold(d,q,fill::randu);
  arma::mat Wnew(d,q,fill::zeros);
  double sig2old = 1.0;
  double sig2new = 1.0;
  arma::vec alpha(q,fill::randu);
  arma::mat Iq(q,q,fill::eye);
  arma::mat Mold = (Wold.t()*Wold) + (sig2old*Iq);
  arma::mat Mnew(q,q,fill::zeros);

  arma::mat MinvWt(q,d,fill::zeros); // precomputing inv(M)*t(W)
  arma::mat emXn(q,N,fill::zeros);
  arma::mat W_LHS(d,q,fill::zeros);  // W : outer products
  arma::mat W_RHS(q,q,fill::zeros);  // W : before pinv

  double term1, term2, term3;

  // 10-3. iterate
  int itercount = 0;
  double tolgap = 100.0;
  arma::vec stopper(2,fill::zeros);
  while (tolgap > reltol){
    // 10-3-1. precompute
    MinvWt = solve(Mold,Wold.t()); // precomputing inv(M)*t(W)
    for (int i=0;i<N;i++){
      emXn.col(i) = MinvWt*T.col(i);
    }

    // 10-3-2. Wnew :: Update W
    W_LHS.zeros(); // clear up previous values
    for (int i=0;i<N;i++){
      W_LHS = W_LHS + auxiliary_outer(T.col(i), emXn.col(i));
    }
    W_RHS = (sig2old*arma::diagmat(alpha)) + sig2old*(static_cast<double>(N))*Mold; // this part N
    for (int i=0;i<N;i++){
      W_RHS = W_RHS + auxiliary_outer(emXn.col(i), emXn.col(i));
    }
    Wnew = W_LHS*pinv(W_RHS);
    Mnew = (Wnew.t()*Wnew) + (sig2old*Iq);

    // 10-3-3. sig2new :: Update Sigma
    sig2new = 0.0;
    for (int i=0;i<N;i++){
      term1 = std::pow(arma::norm(T.col(i),2),2);
      term2 = -2*arma::dot(emXn.col(i), Wnew.t()*T.col(i));
      term3 = arma::trace((sig2old*Mold + auxiliary_outer(emXn.col(i), emXn.col(i)))*Wnew.t()*Wnew);
      sig2new = sig2new + (term1+term2+term3);
    }
    sig2new = sig2new/(static_cast<double>(d*N));

    // 10-3-4. update alpha
    for (int i=0;i<q;i++){
      alpha(i) = static_cast<double>(d)/std::pow(arma::norm(Wnew.col(i)), 2);
    }

    // 10-3-5. stopping criterion
    itercount += 1;
    stopper(0) = norm(Wnew-Wold);
    stopper(1) = std::sqrt(std::pow(sig2new-sig2old,2));
    tolgap = arma::min(stopper);

    // 10-3-6. update
    Wold = Wnew;
    Mold = Mnew;
    sig2old = sig2new;
    if (itercount > maxiter){
      break;
    }
  }

  // 10-4. return list of objects
  return Rcpp::List::create(Rcpp::Named("W")=Wold,
                            Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("sig2")=sig2old,
                            Rcpp::Named("itercount")=itercount
                            );
}

/*
 * 11. Extended Locality Preserving Projection
 * NOTE : simply Z-transform via the S_ij mapping
 */
//' @keywords internal
// [[Rcpp::export]]
arma::mat method_trfextlpp(arma::mat& D, double a, double b){
  const int n = D.n_rows;
  arma::mat output(n,n,fill::zeros);
  double tval = 0.0;
  double outval = 0.0;
  for (int i=0;i<(n-1);i++){
    for (int j=i;j<n;j++){
      tval = D(i,j);

      if (tval<=a){
        outval = 1.0;
      } else if ((a<tval)&&(tval<=((a+b)/2.0))){
        outval = 1-2*std::pow((tval-a)/(b-a), 2.0);
      } else if ((((a+b)/2)<tval)&&(tval<=b)){
        outval = 2*std::pow((tval-a)/(b-a),2);
      } else {
        outval = 0.0;
      }

      output(i,j) = outval;
      output(j,i) = outval;
    }
  }
  return(output);
}

/*
 * 12. Local Similarity Preserving Projection
 */
//' @keywords internal
// [[Rcpp::export]]
arma::mat method_lspp_computeW(arma::mat& S, arma::vec& svec){
  const int n = S.n_rows;
  arma::mat W(n,n,fill::zeros);
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      if (S(i,j)>=svec(i)){
        W(i,j) = S(i,j);
      }
    }
  }
  return(W);
}


/*
 * 13. Kernelized Maximum Margin Criterion
 */
//' @keywords internal
// [[Rcpp::export]]
arma::vec method_kmmcvec(arma::mat& X, arma::mat& partmat, double param){
  // 1. get parameter
  const int n  = X.n_rows;
  const int ni = partmat.n_rows;
  double param22 = 2*param*param;
  // 2. compute
  arma::rowvec vec1(n,fill::zeros);
  arma::rowvec vec2(n,fill::zeros);
  arma::vec output(n,fill::zeros);
  for (int i=0;i<n;i++){
    double tmp = 0.0;
    vec1 = X.row(i);
    for (int j=0;j<ni;j++){
      vec2 = partmat.row(j);
      tmp += exp(-(norm(vec1-vec2,2)*norm(vec1-vec2,2))/param22);
    }
    output(i) = tmp/ni;
  }
  return(output);
}

/*
 * 14. Local Fisher Discriminant Analysis
 */
//' @keywords internal
// [[Rcpp::export]]
double method_lfda_maximaldistance(arma::rowvec& tvec, arma::mat& tmat){
  // 1. get parameter
  const int N = tmat.n_rows;
  // const int p = tvec.n_elem; // unused flag

  // 2. output : iteratively update
  double output = 0.0;
  double tmpout = 0.0;
  // 3. variable declaration 64
  arma::rowvec vec1;
  arma::rowvec vec2;
  vec1 = tvec;
  for (int i=0;i<N;i++){
    vec2 = tmat.row(i);
    tmpout = norm(vec1-vec2,2);
    if (tmpout > output){
      output = tmpout;
    }
  }
  return(output);
}


/*
 * 15. NNPROJMAX & NNPROJMIN
 *            Used in Nonnegative Projection Methods in LINEAR.49.
 *            Note that this function will be utilized further.
 */
//' @keywords internal
// [[Rcpp::export]]
arma::mat method_nnprojmax(arma::mat& C, arma::mat& Uinit, const double tol, const int maxiter){
  // 1. size information
  const int p = Uinit.n_rows;
  const int ndim = Uinit.n_cols;

  // 2. initialize variables;
  arma::mat Uold = Uinit;
  arma::mat Unew(p,ndim,fill::zeros);

  arma::mat C1(p,p,fill::zeros); // for C+
  arma::mat C2(p,p,fill::zeros); // for C-
  C1 = handy_plus(C);
  C2 = C1-C;

  arma::mat term1(p,ndim,fill::zeros); // for numerator   term
  arma::mat term2(p,ndim,fill::zeros); // for denominator term

  // 3. iterate
  double incstop = 100.0;
  int iter = 0;
  while (incstop > tol){
    // 3-1. numerator term
    term1 = C1*Uold + Uold*Uold.t()*C2*Uold;
    term2 = C2*Uold + Uold*Uold.t()*C1*Uold;

    // 3-2. update U
    Unew = handy_hadamardABCsqrt(Uold, term1, term2);

    // 3-3. update incstop
    incstop = arma::norm(Uold-Unew,"fro")/arma::norm(Uold,"fro");

    // 3-4. update U matrix itself
    Uold = Unew;

    // 3-5. iteration numbers
    iter = iter + 1;
    if (iter >= maxiter){
      break;
    }
  }

  // 4. return output
  return(Uold);
}
//' @keywords internal
// [[Rcpp::export]]
arma::mat method_nnprojmin(arma::mat& C, arma::mat& Uinit, const double tol, const int maxiter){
  // 1. size information
  const int p = Uinit.n_rows;
  const int ndim = Uinit.n_cols;

  // 2. initialize variables;
  arma::mat Uold = Uinit;
  arma::mat Unew(p,ndim,fill::zeros);

  arma::mat C1(p,p,fill::zeros); // for C+
  arma::mat C2(p,p,fill::zeros); // for C-
  C1 = handy_plus(C);
  C2 = C1-C;

  arma::mat term1(p,ndim,fill::zeros); // for numerator   term
  arma::mat term2(p,ndim,fill::zeros); // for denominator term

  // 3. iterate
  double incstop = 100.0;
  int iter = 0;
  while (incstop > tol){
    // 3-1. numerator term
    term1 = C2*Uold + Uold*Uold.t()*C1*Uold;
    term2 = C1*Uold + Uold*Uold.t()*C2*Uold;

    // 3-2. update U
    Unew = handy_hadamardABCsqrt(Uold, term1, term2);

    // 3-3. update incstop
    incstop = arma::norm(Uold-Unew,"fro")/arma::norm(Uold,"fro");

    // 3-4. update U matrix itself
    Uold = Unew;

    // 3-5. iteration numbers
    iter = iter + 1;
    if (iter >= maxiter){
      break;
    }
  }

  // 4. return output
  return(Uold);
}



/*
 * 16. NNEMBEDMIN
 *            Used in solving tr(YMY') s.t YY'=I, Y>=0
 */
//' @keywords internal
// [[Rcpp::export]]
arma::mat method_nnembedmin(arma::mat& M, arma::mat& Yinit, const double tol, const int maxiter){
  // 1. size information
  const int m = Yinit.n_rows;
  const int n = Yinit.n_cols;

  // 2. initialize variables;
  arma::mat Yold = Yinit;
  arma::mat Ynew(m,n,fill::zeros);

  arma::mat M1(n,n,fill::zeros); // for M+
  arma::mat M2(n,n,fill::zeros); // for M-
  M1 = handy_plus(M);
  M2 = M1-M;

  arma::mat term1(m,n,fill::zeros); // for numerator   term
  arma::mat term2(m,n,fill::zeros); // for denominator term

  // 3. iterate
  double incstop = 100.0;
  int iter = 0;
  while (incstop > tol){
    // 3-1. numerator term
    term1 = Yold*M2 + Yold*M1*Yold.t()*Yold;
    term2 = Yold*M1 + Yold*M2*Yold.t()*Yold;

    // 3-2. update U
    Ynew = handy_hadamardABCsqrt(Yold, term1, term2);

    // 3-3. update incstop
    incstop = arma::norm(Yold-Ynew,"fro")/arma::norm(Yold,"fro");

    // 3-4. update U matrix itself
    Yold = Ynew;

    // 3-5. iteration numbers
    iter = iter + 1;
    if (iter >= maxiter){
      break;
    }
  }

  // 4. return output
  return(Yold);
}

/*
 * 17. SPUFS : return score of weight vectors
 */
//' @keywords internal
// [[Rcpp::export]]
arma::vec method_spufs(arma::mat& X, arma::mat Ls, double alpha, double beta, double epsilon){
  // basic setup
  const int n = X.n_rows;
  const int d = X.n_cols;

  // initialize
  arma::mat Qold(d,d,fill::eye);
  arma::mat Qnew(d,d,fill::zeros);

  arma::mat Wold(d,d,fill::eye);
  arma::mat Wnew(d,d,fill::zeros);

  arma::mat XtX   = X.t()*X;
  arma::mat LHS1  = (beta*X.t()*Ls*X + XtX); // within LHS, it's fixed.
  arma::mat RHS   = XtX;

  arma::mat LHS(d,d,fill::zeros);
  double increment = 100000.0;

  // main iteration : stopping criterion Q
  arma::colvec wcoli;
  while (increment > 1e-5){
    // 1. update W
    LHS  = LHS1 + alpha*Qold;
    Wnew = arma::solve(LHS, RHS);

    // 2. update Q
    for (int i=0;i<d;i++){
      wcoli     = Wnew.col(i);
      Qnew(i,i) = 1.0/(2.0*std::sqrt(arma::as_scalar(wcoli.t()*wcoli) + epsilon));
    }

    // 3. update increment
    increment = arma::norm(Qold-Qnew);

    // 4. update Qold and Qnew
    Wold = Wnew;
    Qold = Qnew;
  }

  // compute score for each
  arma::vec output(d,fill::zeros);
  for (int i=0;i<d;i++){
    wcoli = Wold.col(i);
    output(i) = arma::norm(wcoli, 2);
  }
  return(output);
}

