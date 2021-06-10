#' Laplacian Eigenmaps
#'
#' \code{do.lapeig} performs Laplacian Eigenmaps (LE) to discover low-dimensional
#' manifold embedded in high-dimensional data space using graph laplacians. This
#' is a classic algorithm employing spectral graph theory.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param weighted \code{TRUE} for weighted graph laplacian and \code{FALSE} for
#' combinatorial laplacian where connectivity is represented as 1 or 0 only.
#' @param kernelscale kernel scale parameter. Default value is 1.0.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{eigvals}{a vector of eigenvalues for laplacian matrix.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' lab   = as.factor(iris[subid,5])
#'
#' ## try different levels of connectivity
#' out1 <- do.lapeig(X, type=c("proportion",0.10), weighted=FALSE)
#' out2 <- do.lapeig(X, type=c("proportion",0.20), weighted=FALSE)
#' out3 <- do.lapeig(X, type=c("proportion",0.50), weighted=FALSE)
#'
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=lab, main="10% connected")
#' plot(out2$Y, pch=19, col=lab, main="20% connected")
#' plot(out3$Y, pch=19, col=lab, main="50% connected")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{belkin_laplacian_2003}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_LAPEIG
#' @concept nonlinear_methods
#' @export
do.lapeig <- function(X, ndim=2, type=c("proportion",0.1),
                      symmetric=c("union","intersect","asymmetric"),
                      preprocess=c("null","center","scale","cscale","whiten","decorrelate"),
                      weighted=FALSE, kernelscale=1.0){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lapeig : 'ndim' is a positive integer in [1,#(covariates)).")}

  # 2. ... parameters
  # 2-1. aux.graphnbd
  #   1. type        : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #   2. symmetric   : 'intersect','union', or 'asymmetric'
  # 2-2. Eigenmaps
  #   1. preprocess  : 'null', 'center','decorrelate', or 'whiten'
  #   2. weighted    : TRUE (default) / FALSE (binary)
  #   3. kernelscale : Kernel Scale Parameter (default; 1.0) / Infty -> binary

  nbdtype = type
  if (missing(symmetric)){
    nbdsymmetric = "union"
  } else {
    nbdsymmetric = match.arg(symmetric)
  }
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  wflag = weighted
  if (!is.logical(wflag)){
    stop("* do.lapeig : 'weighted' flag should be a logical input.")
  }
  t = kernelscale
  if (!is.numeric(t)||is.na(t)||(t<=0)){
    stop("* do.lapeig : 'kernelscale' is a positive real value.")
  }
  if (t==Inf){
    wflag  = FALSE
  }

  # 3. process : data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  n = nrow(pX)
  p = ncol(pX)

  # 4. process : neighborhood selection
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)

  # 5. main computation
  #   5-1. Compute Weight Matrix W
  if (wflag==FALSE){
    W = matrix(as.double(nbdstruct$mask),nrow=nrow(nbdstruct$mask))
    W = (W+t(W))/2
  } else {
    W = exp((-(nbdstruct$mask*nbdstruct$dist)^2)/t)
    idxnan = which(is.nan(W))
    W[idxnan] = 0
    diag(W) = 0
    W = (W+t(W))/2
  }
  #   5-2. Compute Embedding
  embedding = method_eigenmaps(W);

  # 6. Output
  #   this uses lowest (ndim+1) eigenpairs
  eigvals = embedding$eigval
  eigvecs = embedding$eigvec

  result = list()
  result$Y = eigvecs[,2:(ndim+1)]
  result$eigvals = eigvals[2:(ndim+1)]
  trfinfo$algtype = "nonlinear"
  result$trfinfo  = trfinfo
  return(result)
}
