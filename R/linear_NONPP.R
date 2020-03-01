#' Nonnegative Orthogonal Neighborhood Preserving Projections
#'
#' Nonnegative Orthogonal Neighborhood Preserving Projections (NONPP) is a variant of ONPP where
#' projection vectors - or, basis for learned subspace - contain no negative values.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center" and other options of "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param maxiter number of maximum iteraions allowed.
#' @param reltol stopping criterion for incremental relative error.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \dontrun{
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])+50
#' label = as.integer(iris$Species)
#'
#' ## use different levels of connectivity
#' out1 = do.nonpp(X, type=c("proportion",0.1))
#' out2 = do.nonpp(X, type=c("proportion",0.2))
#' out3 = do.nonpp(X, type=c("proportion",0.5))
#'
#' ## visualize
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(out1$Y, col=label, main="NONPP::10% connected")
#' plot(out2$Y, col=label, main="NONPP::20% connected")
#' plot(out3$Y, col=label, main="NONPP::50% connected")
#' par(opar)
#' }
#'
#' @seealso \code{\link{do.onpp}}
#' @references
#' \insertRef{zafeiriou_nonnegative_2010}{Rdimtools}
#'
#' @rdname linear_NONPP
#' @author Kisung You
#' @export
do.nonpp <- function(X, ndim=2, type=c("proportion",0.1),
                     preprocess=c("null","center","decorrelate","whiten"),
                     maxiter=1000, reltol=1e-5){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.nonpp : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   3. type
  nbdtype = type
  nbdsymmetric = "union"
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   * maxiter and reltol
  maxiter = as.integer(maxiter)
  if (!check_NumMM(maxiter, 5, 1e+6)){stop("* do.nonpp : 'maxiter' is a large positive integer for the number of iterations.")}
  reltol = as.double(reltol)
  if (!check_NumMM(reltol, .Machine$double.eps, 1)){stop("* do.nonpp : 'reltol' is a small positive real number for stopping criterion.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY and LLE step
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)

  #   3. LLE computation
  regparam = 1.0
  W = array(0,c(n,n))
  for (i in 1:n){
    #   3-1. separate target mask vector
    tgtidx  = which(nbdstruct$mask[i,])
    #   3-2. select data
    #        For convenience, target matrix is transposed for Armadillo
    vec_tgt = pX[i,]
    mat_tgt = t(pX[tgtidx,])
    k = ncol(mat_tgt)
    #   3-3. no automatic regularization
    W[i,tgtidx] = method_lleW(mat_tgt,vec_tgt,regparam);
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION FOR N-ONPP
  #   1. cost function : C
  diagN  = diag(n)
  C      = t(pX)%*%(t(diagN-W)%*%(diagN-W))%*%pX
  #   2. initialize for U
  Uinit = matrix(runif(p*ndim),nrow=p)
  #   3. compute projection matrix
  projection = method_nnprojmin(C, Uinit, reltol, maxiter)
  #   4. additional step : NA
  projection[(is.na(projection)||(is.infinite(projection)))] = 1
  for (i in 1:ndim){
    tgt = as.vector(projection[,i])
    projection[,i] = tgt/sqrt(sum(tgt^2))
  }
  projection = aux.adjprojection(projection)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
