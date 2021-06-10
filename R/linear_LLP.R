#' Local Learning Projections
#'
#' While Principal Component Analysis (PCA) aims at minimizing global estimation error, Local Learning
#' Projection (LLP) approach tries to find the projection with the minimal \emph{local}
#' estimation error in the sense that each projected datum can be well represented
#' based on ones neighbors. For the kernel part, we only enabled to use
#' a gaussian kernel as suggested from the original paper. The parameter \code{lambda}
#' controls possible rank-deficiency of kernel matrix.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}.
#' See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param t bandwidth for heat kernel in \eqn{(0,\infty)}.
#' @param lambda regularization parameter for kernel matrix in \eqn{[0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \donttest{
#' ## generate data
#' set.seed(100)
#' X <- aux.gensamples(n=100, dname="crown")
#'
#' ## test different lambda - regularization - values
#' out1 <- do.llp(X,ndim=2,lambda=0.1)
#' out2 <- do.llp(X,ndim=2,lambda=1)
#' out3 <- do.llp(X,ndim=2,lambda=10)
#'
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, main="lambda=0.1")
#' plot(out2$Y, pch=19, main="lambda=1")
#' plot(out3$Y, pch=19, main="lambda=10")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{wu_local_2007}{Rdimtools}
#'
#' @rdname linear_LLP
#' @concept linear_methods
#' @export
do.llp <- function(X, ndim=2, type=c("proportion",0.1), symmetric=c("union","intersect","asymmetric"),
                   preprocess = c("center","scale","cscale","decorrelate","whiten"), t=1.0, lambda=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.llp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. type
  nbdtype = type
  #   4. symmetric
  nbdsymmetric = match.arg(symmetric)
  #   5. preprocess
  algpreprocess = match.arg(preprocess)
  #   6. t : kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t,0,1e+10,compact=FALSE)){stop("* do.llp : 't' should be a positive real number.")}
  #   7. lambda : adjusting kernel matrix
  lambda = as.double(lambda)
  if (!check_NumMM(lambda,0,1e+10,compact=TRUE)){stop("* do.llp : 'lambda' should be a nonnegative real number.")}

  #------------------------------------------------------------------------
  ## COMPUTATION PART 1 : PREPROCESSING
  #   1. preprocessing the data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. find neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  pXmask    = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION PART 2 : LLP
  #   Compute A
  A = array(0,c(n,n))
  for (i in 1:n){
    #   1. target vector X_i
    tgtvec = pX[i,]
    #   2. target matrix X_j's
    tgtidx = setdiff(which(pXmask[i,]),i)
    tgtmat = pX[tgtidx,]
    #   3. compute ki as row-1 matrix (rowvec)
    tmp_ki = llp_compute_ki(tgtvec, tgtmat, t)
    #   4. compute Ki
    tmp_Ki = llp_compute_Ki(tgtmat, t, lambda)
    #   5. compute and assign values for A
    A[i,tgtidx] = as.vector(tmp_ki%*%(base::solve(tmp_Ki)))
  }
  #   Compute T
  IA   = (diag(n)-A)
  matT = (t(IA)%*%IA)
  #   Projection
  XTX  = t(pX)%*%matT%*%pX

  #------------------------------------------------------------------------
  ## RETURN OUTPUT
  projection = aux.adjprojection(base::eigen(XTX)$vectors[,p:(p-ndim+1)])

  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}








#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
llp_compute_ki <- function(xi,matxj,t){
  n = nrow(matxj)
  vecki = rep(0,n)
  for (i in 1:n){
    tgtvec = as.vector(xi)-as.vector(matxj[i,])
    vecki[i] = exp(-sum(vecki*vecki)/(2*(t^2)))
  }
  vecki = matrix(vecki,nrow=1)
  return(vecki)
}
#' @keywords internal
#' @noRd
llp_compute_Ki <- function(matxj, t, lambda){
  n = nrow(matxj)
  dmat = exp(-(as.matrix(dist(matxj))^2)/(2*(t^2)))
  Ki = dmat + lambda*diag(n)
  return(Ki)
}
