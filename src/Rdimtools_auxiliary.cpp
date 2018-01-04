/*
 * 1. aux_preprocess : center, decorrelate, or whiten
 * 2. aux_perplexity : given target perplexity, compute P
 * 3. aux_shortestpath : Floyd-Warshall algorithm.
 * 4. aux_landmarkMaxMin : select landmark points using MaxMin tactic
 * 5. aux_kernelcov : compute K and centered K matrix
 * 6. aux_eigendecomposition : eigendecomposition of a given symmetric matrix
 * 7. aux_minmax : find minimum and maximum values for each dimension
 * 8. aux_regout : regress out a vector on a matrix : row-sense
 * 9. aux_scatter          : sum{(x_i-mu)(x_i-mu)^T}
 *    aux_scatter_pairwise : sum{sum{(x_i-x_j)(x_i-x_j)^T}}
 *
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// subroutine1
double computeH(arma::vec& D, const double var, const int i, const int n){
  vec pji(n);
  for (int j=0;j<n;j++){
    if (j!=i){
      pji(j) = exp(-D(j)/(2*var));
    } else{
      pji(j) = 0;
    }
  }
  double cdenom = 0;
  for (int j=0;j<n;j++){
    if (j!=i){
      cdenom += exp(-D(j)/(2*var));
    }
  }
  pji /= cdenom;

  double H = 0;
  for (int j=0;j<n;j++){
    H -= pji(j)*log2(pji(j));
  }
  return(H);
}

// subroutine2
LogicalMatrix isweird(NumericMatrix x){
  const int n = x.nrow();
  LogicalMatrix out(n,n);

  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      out(i,j) = ((x(i,j)==R_NegInf)||(x(i,j)==R_PosInf)||(NumericVector::is_na(x(i,j))));
    }
  }
  return out;
}

/*
 * 1. aux_preprocess : center, decorrelate, or whiten
 */
// [[Rcpp::export]]
Rcpp::List aux_preprocess(arma::mat& X, const int flag){
  const int p = X.n_rows;
  const int n = X.n_cols;
  arma::vec Xmean = mean(X,1);

  // 1-1. centering
  for (int i=0;i<n;i++){
    X.col(i) -= Xmean;
  }
  if (flag==1){
    return Rcpp::List::create(Rcpp::Named("type")="center",
                              Rcpp::Named("output")=X,
                              Rcpp::Named("mean")=Xmean,
                              Rcpp::Named("multiplier")=1);
  } else if (flag==2){
  // 1-2. decorrelate
    mat covX = (X*(X.t()))/(n-1);
    vec eigval;
    mat eigvec;
    eig_sym(eigval,eigvec,covX);

    mat output = (eigvec.t())*X;
    return Rcpp::List::create(Rcpp::Named("type")="decorrelate",
                              Rcpp::Named("output")=output,
                              Rcpp::Named("mean")=Xmean,
                              Rcpp::Named("multiplier")=eigvec);
  } else if (flag==3){
  // 1-3. whiten
    mat covX = (X*X.t())/(n-1);
    vec eigval;
    mat eigvec;
    eig_sym(eigval,eigvec,covX);

    mat weight = diagmat(1/sqrt(eigval));
    mat multiplier = eigvec*(weight.t());
    mat output = multiplier.t()*X;
    return Rcpp::List::create(Rcpp::Named("type")="whiten",
                              Rcpp::Named("output")=output,
                              Rcpp::Named("mean")=Xmean,
                              Rcpp::Named("multiplier")=multiplier);
  } else {
    Rcpp::stop("choose either one of three options.");
    return 0;
  }
}

/*
 * 2. aux_perplexity : given target perplexity, compute P
 */
// [[Rcpp::export]]
Rcpp::List aux_perplexity(arma::mat& X,const double perplexity){
  // 2-1. basic settings
  const int d = X.n_rows;
  const int n = X.n_cols;
  const double tol = 1e-5;
  const double logU= log2(perplexity);

  // 2-2. initialize
  mat P = zeros<mat>(n,n);
  vec beta = ones<vec>(n);
  mat D(n,n);
  for (int i=0;i<n;i++){
    for (int j=(i+1);j<n;j++){
      double dval = std::pow(norm(X.col(i)-X.col(j)),2);
      D(i,j) = dval;
      D(j,i) = dval;
    }
  }

  // 2-3. run over all datapoints
  for (int i=0;i<n;i++){
    // 2-3-1. set minimum and maximum for precision
    double betamin = -1e+10;
    double betamax = 1e+10;

    // 2-3-2. initialize
    vec Dvec = D.col(i);
    double H = computeH(Dvec,beta(i),i,n);

    double Hdiff = H - logU;
    int tries = 0;
    while ((std::abs(Hdiff) > tol)&&(tries<50)){
      if (Hdiff > 0){
        betamin = beta(i);
        beta(i) = (beta(i)+betamax)/2;
      } else{
        betamax = beta(i);
        beta(i) = (beta(i)+betamin)/2;
      }
      // 2-3-3. recompute values
      H = computeH(Dvec,beta(i),i,n);
      Hdiff = H - logU;
      tries += 1;
    }
  }

  // 2-4. now we have all the values - just compute it again.
  for (int i=0;i<n;i++){
    double cdenom = 0;
    double pval2 = 2*std::pow(beta(i),2);
    for (int k=0;k<n;k++){
      if (k!=i){
        cdenom += exp(-D(i,k)/pval2);
      }
    }
    for (int j=0;j<n;j++){
      if (j!=i){
        P(i,j) = (exp(-D(i,j)/pval2))/cdenom;
      }
    }
  }

  // 2-5. return output
  return Rcpp::List::create(Rcpp::Named("vars")=beta,
                     Rcpp::Named("P")=P);
}

/*
 * 3. aux_shortestpath : Floyd-Warshall algorithm.
 */
// [[Rcpp::export]]
Rcpp::NumericMatrix aux_shortestpath(NumericMatrix& wmat){
  // 3-1. get ready
  const int v = wmat.nrow();
  NumericMatrix dist(v,v);
  for (int i=0;i<v;i++){
    for (int j=0;j<v;j++){
      dist(i,j) = R_PosInf;
    }
  }
  // 3-2. initialization
  LogicalMatrix checker = isweird(wmat);

  // 3-3. Floyd-Warshall algorithm
  // 3-3-1. vertex
  for (int i=0;i<v;i++){
    dist(i,i) = 0;
  }
  // 3-3-2. edge list
  for (int i=0;i<v;i++){
    for (int j=0;j<v;j++){
      if (checker(i,j)==false){
        dist(i,j) = wmat(i,j);
      }
    }
  }
  // 3-3-3. main iteration
  for (int k=0;k<v;k++){
    for (int i=0;i<v;i++){
      for (int j=0;j<v;j++){
        if (dist(i,j)>(dist(i,k)+dist(k,j))){
          dist(i,j)=dist(i,k)+dist(k,j);
        }
      }
    }
  }

  // 3-4. return output
  return(dist);
}

/*
 * 4. aux_landmarkMaxMin : select landmark points using MaxMin tactic
 *    Note that those vectors of indices should be adjusted as well as the result
 */
// [[Rcpp::export]]
int aux_landmarkMaxMin(arma::mat& pD, arma::vec& plandmark, arma::vec& seqnp){
  // 4-1. basic setting
  const int nlandmark = plandmark.n_elem;
  const int ntestpts  = seqnp.n_elem;

  // 4-2. we should be careful ; -1 for both vectors
  vec veclandmark = plandmark - 1;
  vec vecseqnp    = seqnp - 1;

  // 4-3. main iteration
  int currentidx      = 0;
  double currentdists = 123456789;
  for (int i=0;i<ntestpts;i++){
    int testpt       = vecseqnp(i);
    double testdists = 0;
    for (int j=0;j<nlandmark;j++){
      int targetpt = veclandmark(j);
      testdists += pD(testpt,targetpt);
    }
    if (testdists<currentdists){
      currentidx   = testpt;
      currentdists = testdists;
    }
  }
  currentidx += 1;

  // 4-4. return output
  return(currentidx);
}

/*
 * 5. aux_kernelcov : compute K and centered K matrix
 */
typedef double (*kernelPtr)(arma::vec& x, arma::vec& y, const double par1, const double par2);
double kernel_linear(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = as_scalar(dot(x,y))+par1;
  return(result);
}
double kernel_polynomial(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = pow(as_scalar(dot(x,y))+par1,par2);
  return(result);
}
double kernel_gaussian(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = exp(pow(norm(x-y),2)*(-par1));
  return(result);
}
double kernel_laplacian(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = exp((-par1)*norm(x-y));
  return(result);
}
double kernel_anova(arma::vec& x, arma::vec& y, const double par1, const double par2){
  const int n = x.n_elem;
  double result = 0;
  for (int i=0;i<n;i++){
    result += pow(exp((-par1)*pow((x(i)-y(i)),2)),par2);
  }
  return(result);
}
double kernel_sigmoid(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = tanh((as_scalar(dot(x,y)))*par1 +par2);
  return(result);
}
double kernel_rq(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double n2 = as_scalar(pow(norm(x-2),2));
  double result = 1-(n2/(n2+par1));
  return(result);
}
double kernel_mq(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = (pow(norm(x-y),2) + pow(par1,2));
  return(result);
}
double kernel_iq(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = 1/(pow(norm(x-y),2)+pow(par1,2));
  return(result);
}
double kernel_imq(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = sqrt(1/(pow(norm(x-y),2)+pow(par1,2)));
  return(result);
}
double kernel_circular(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double nxy = norm(x-y);
  double result;
  if (nxy < par1){
    result = (2/par2)*acos(-nxy/par1) - (2/par2)*(nxy/par1)*sqrt(1-pow(nxy/par1,2));
  } else {
    result = 0;
  }
  return(result);
}
double kernel_spherical(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double nxy = norm(x-y);
  double result;
  if (nxy < par1){
    result = 1-(1.5*nxy/par1)+(0.5*pow(nxy/par1,3));
  } else {
    result = 0;
  }
  return(result);
}
double kernel_power(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = -pow(norm(x-y),par1);
  return(result);
}
double kernel_log(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = -log(pow(norm(x-y),par1)+1);
  return(result);
}
double kernel_spline(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result=1;
  const int n = x.n_elem;
  double tx,ty,mxy;
  for (int i=0;i<n;i++){
    tx = x(i);
    ty = y(i);
    if (tx<ty){
      mxy = tx;
    } else {
      mxy = ty;
    }
    result *= 1+(tx*ty)+(tx*ty*mxy)-(((tx+ty)/2)*pow(mxy,2))+(pow(mxy,3)/3);
  }
  return(result);
}
double kernel_cauchy(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = pow(par1,2)/(pow(par1,2) + pow(norm(x-y),2));
  return(result);
}
double kernel_chisq(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = 0;
  const int n = x.n_elem;
  for (int i=0;i<n;i++){
    result += (x(i)*y(i)*2)/(x(i)+y(i));
  }
  return(result);
}
double kernel_histintx(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = 0;
  const int n = x.n_elem;
  for (int i=0;i<n;i++){
    if (x(i)<y(i)){
      result += x(i);
    } else {
      result += y(i);
    }
  }
  return(result);
}
double kernel_ghistintx(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = 0;
  const int n = x.n_elem;
  double xa,yb,xaa,ybb;
  for (int i=0;i<n;i++){
    xa = x(i);
    yb = y(i);
    if (xa < 0){
      xa *= -1;
    }
    if (yb < 0){
      yb *= -1;
    }
    xaa = pow(xa,par1);
    ybb = pow(yb,par2);
    if (xaa < ybb){
      result += xaa;
    } else {
      result += ybb;
    }
  }
  return(result);
}
double kernel_t(arma::vec& x, arma::vec& y, const double par1, const double par2){
  double result = 1/(1+pow(norm(x-y),par1));
  return(result);
}

XPtr<kernelPtr> decidePtr(const int n){
  if (n==1){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_linear)));
  } else if (n==2){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_polynomial)));
  } else if (n==3){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_gaussian)));
  } else if (n==4){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_laplacian)));
  } else if (n==5){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_anova)));
  } else if (n==6){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_sigmoid)));
  } else if (n==7){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_rq)));
  } else if (n==8){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_mq)));
  } else if (n==9){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_iq)));
  } else if (n==10){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_imq)));
  } else if (n==11){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_circular)));
  } else if (n==12){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_spherical)));
  } else if (n==13){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_power)));
  } else if (n==14){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_log)));
  } else if (n==15){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_spline)));
  } else if (n==16){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_cauchy)));
  } else if (n==17){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_chisq)));
  } else if (n==18){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_histintx)));
  } else if (n==19){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_ghistintx)));
  } else if (n==20){
    return(XPtr<kernelPtr>(new kernelPtr(&kernel_t)));
  } else {
    return XPtr<kernelPtr>(R_NilValue);
  }
}
// [[Rcpp::export]]
Rcpp::List aux_kernelcov(arma::mat& tX, const int knumber, const double par1, const double par2){
  // 5-1. basic settings
  const int N = tX.n_cols;
  mat K(N,N,fill::zeros);
  mat onesN(N,N,fill::ones);
  mat Kcenter(N,N,fill::zeros);

  // 5-2. get pointer
  XPtr<kernelPtr> xpfun = decidePtr(knumber);
  kernelPtr trfK = *xpfun;

  // 5-3. main iteration
  vec row1;
  vec row2;
  double tmpcompute;
  for (int i=0;i<N;i++){
    row1 = tX.col(i);
    for (int j=i;j<N;j++){
      if (i==j){
        K(i,i) = as_scalar(trfK(row1,row1,par1,par2));
      } else {
        row2 = tX.col(j);
        tmpcompute = trfK(row1,row2,par1,par2);
        K(i,j) = tmpcompute;
        K(j,i) = tmpcompute;
      }
    }
  }

  // 5-4. centering of K
  double constN2 = std::pow(N,2);
  Kcenter = K - ((onesN*K)/N) - ((K*onesN)/N) + (onesN*K*onesN)/constN2;

  // 5-5. return output
  return Rcpp::List::create(Rcpp::Named("K")=K,
                            Rcpp::Named("Kcenter")=Kcenter);
}


/*
 * 6. aux_eigendecomposition : eigendecomposition of a given symmetric matrix
 */
// [[Rcpp::export]]
Rcpp::List aux_eigendecomposition(arma::mat& X){
  arma::vec eigval;
  arma::mat eigvec;
  arma::mat sX = ((X+X.t())/2);

  arma::eig_sym(eigval, eigvec, sX);
  return Rcpp::List::create(Rcpp::Named("eigval")=eigval,
                            Rcpp::Named("eigvec")=eigvec);
}

/*
 * 7. aux_minmax : find minimum and maximum values for each dimension
 */
// [[Rcpp::export]]
arma::mat aux_minmax(arma::mat& X, const double gap){
  // 7-1. basic setting
  const int n = X.n_rows;
  const int d = X.n_cols;
  mat output(2,d,fill::zeros);

  // 7-2. iteration
  for (int i=0;i<d;i++){
    output(0,i) = as_scalar(X.col(i).min()) - gap;
    output(1,i) = as_scalar(X.col(i).max()) + gap;
  }

  // 7-3. return output
  return(output);
}

/*
 * 8. aux_regout : regress out a vector on a matrix : row-sense
 */
// [[Rcpp::export]]
arma::mat aux_regout(arma::mat& X, arma::rowvec tgt){
  // 8-1. basic setting
  const int n = X.n_rows;
  const int p = X.n_cols;
  arma::mat output(n,p,fill::zeros);

  // 8-2. main iteration
  // Note that vector 'tgt' should be a length of p.
  for (int i=0;i<n;i++){
    output.row(i) = X.row(i) - (dot(X.row(i),tgt))*tgt;
  }

  // 8-3. return output
  return(output);
}

/*
 * 9-1. aux_scatter          : sum{(x_i-mu)(x_i-mu)^T}
 * 9-2. aux_scatter_pairwise : sum{sum{(x_i-x_j)(x_i-x_j)^T}}
 */
// [[Rcpp::export]]
arma::mat aux_scatter(arma::mat& X, arma::rowvec& mu){
  // 1. parameters
  const int n = X.n_rows;
  const int p = X.n_cols;
  // 2. output and settings
  arma::mat output(p,p,fill::zeros);

  arma::rowvec vecrow(p,fill::zeros);
  arma::colvec veccol(p,fill::zeros);
  // 3. iteration
  for (int i=0;i<n;i++){
    vecrow = X.row(i)-mu;
    veccol = vecrow.t();
    output += (veccol*vecrow);
  }
  return(output);
}
// [[Rcpp::export]]
arma::mat aux_scatter_pairwise(arma::mat& X){
  // 1. parameters
  const int n = X.n_rows;
  const int p = X.n_cols;
  // 2. output and setting
  arma::mat output(p,p,fill::zeros);
  arma::rowvec vec1(p,fill::zeros);
  arma::rowvec vec2(p,fill::zeros);

  arma::colvec veccol(p,fill::zeros);
  arma::rowvec vecrow(p,fill::zeros);
  // 3. iteration
  for (int i=0;i<n;i++){
    vec1 = X.row(i);
    for (int j=0;j<n;j++){
      vec2 = X.row(j);
      if (i!=j){
        vecrow = vec1-vec2;
        veccol = vecrow.t();
        output += (veccol*vecrow);
      }
    }
  }
  return(output);
}

