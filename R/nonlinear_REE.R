#' Robust Euclidean Embedding
#'
#' Robust Euclidean Embedding (REE) is an embedding procedure exploiting
#' robustness of \eqn{\ell_1} cost function. In our implementation, we adopted
#' a generalized version with weight matrix to be applied as well. Its original
#' paper introduced a subgradient algorithm to overcome memory-intensive nature of
#' original semidefinite programming formulation.
#'
#' @param X an \code{(n-by-p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param W an \code{(n-by-n)} weight matrix. Default is uniform weight of 1s.
#' @param preprocess an additional option for preprocessing the data.
#' Default is ``null'' and  other methods of ``center'',``decorrelate'', or ``whiten''
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param initc initial \code{c} value for subgradient iterating stepsize, \eqn{c/\sqrt{i}}.
#' @param dmethod a type of distance measure. See \code{\link[stats]{dist}} for more details.
#' @param maxiter maximum number of iterations for subgradient descent method.
#' @param abstol stopping criterion for subgradient descent method.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \code{(n-by-ndim)} matrix whose rows are embedded observations.}
#' \item{niter}{the number of iterations taken til convergence. }
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#'
#' @examples
#' ## generate swiss roll data
#' X = aux.gensamples(n=123)
#'
#' ## 1. no preprocessing
#' output1 <- do.ree(X,ndim=2,maxiter=50)
#'
#' ## 2. use decorrelated data
#' output2 <- do.ree(X,ndim=2,preprocess="decorrelate",maxiter=50)
#'
#' ## Visualize three different projections
#' par(mfrow=c(1,2))
#' plot(output1$Y[,1],output1$Y[,2],main="centered")
#' plot(output2$Y[,1],output2$Y[,2],main="decorrelated")
#'
#' @references
#' \insertRef{cayton_robust_2006}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_REE
#' @export
do.ree <- function(X, ndim=2, W=NA,preprocess="null",initc=1.0,
                   dmethod=c("euclidean","maximum","manhattan",
                             "canberra","binary","minkowski"),
                   maxiter=100, abstol=1e-3){
  ## PREPROCESSING
  # 1. check X
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  # 2. check ndim
  if (!check_ndim(ndim,p)){
    stop("* do.ree : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  # 3. check W
  if ((missing(W))||(is.na(W))){
    W = matrix(rep(1,n^2),nrow=n)
  } else {
    if (!check_WMAT(W,n)){
      stop("* do.ree : 'W' should be a square matrix of nonnegative real numbers.")
    }
  }
  W = W/sum(W)
  # 4. run preprocessing
  if (preprocess=="null"){
    trfinfo = list()
    trfinfo$type = "null"
    pX = as.matrix(X,nrow=nrow(X))
  } else {
    tmplist = aux.preprocess(X,type=preprocess)
    trfinfo = tmplist$info
    pX      = tmplist$pX
  }
  trfinfo$algtype = "nonlinear"
  # 5. dmethod
  dmethod = match.arg(dmethod)
  # 6. initc
  if ((length(initc)!=1)||(initc<=0)||(is.na(initc))||(is.infinite(initc))){
    stop("* do.ree : 'initc' should be a positive real number.")
  } else {
    initc = 1.0
  }
  # 7. maxiter
  if ((length(maxiter)!=1)||(maxiter<2)||(is.na(maxiter))||(is.infinite(maxiter))||(abs(maxiter-round(maxiter))>sqrt(.Machine$double.eps))){
    stop("* do.ree : 'maxiter' should be a positive integer larger than 1.")
  }
  maxiter = as.integer(maxiter)
  # 8. abstol
  if ((length(abstol)!=1)||(abstol<=sqrt(.Machine$double.eps))||(is.na(abstol))||(is.infinite(abstol))){
    stop("* do.ree : 'abstol' should be a small positive real number.")
  }
  abstol = as.double(abstol)

  ## Main Computation : let's write everything in C
  D = as.matrix(dist(pX, method=dmethod))
  B = matrix(rnorm(n*n),nrow=n); B = t(B)%*%B; # make as SPD
  tmpout = method_ree(B, W, D, initc, abstol, maxiter); # returns $B matrix and $iter the number of iterations

  ## Post Processing for injection : on largest ones
  eigB = eigen(tmpout$B)
  tgtvecs = eigB$vectors[,1:ndim]
  tgteigs = eigB$values[1:ndim]
  tgteigs = ifelse(tgteigs>0,tgteigs,0)

  ## RETURN
  output = list()
  output$Y       = tgtvecs%*%diag(tgteigs)
  output$niter   = tmpout$iter
  output$trfinfo = trfinfo
  return(output)
}
