#' Laplacian Eigenmaps
#'
#' \code{do.lapeig} performs Laplacian Eigenmaps (LE) to discover low-dimensional
#' manifold embedded in high-dimensional data space using graph laplacians. This
#' is a classic algorithm employing spectral graph theory.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param ... extra parameters including \describe{
#' \item{kernelscale}{kernel scale parameter. Default value is 1.0.}
#' \item{preprocess}{an additional option for preprocessing the data.
#' Default is \code{"null"}. See also \code{\link{aux.preprocess}} for more details.}
#' \item{symmetric}{one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.}
#' \item{type}{a vector of neighborhood graph construction. Following types are supported;
#' \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#' Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#' among all data points. See also \code{\link{aux.graphnbd}} for more details.}
#' \item{weighted}{a logical; \code{TRUE} for weighted graph laplacian and \code{FALSE} for
#' combinatorial laplacian where connectivity is represented as 1 or 0 only.}
#' }
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{eigvals}{a vector of eigenvalues for laplacian matrix.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{algorithm}{name of the algorithm.}
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
#' out1 <- do.lapeig(X, type=c("proportion",0.5), weighted=FALSE)
#' out2 <- do.lapeig(X, type=c("proportion",0.10), weighted=FALSE)
#' out3 <- do.lapeig(X, type=c("proportion",0.25), weighted=FALSE)
#'
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=lab, main="5% connected")
#' plot(out2$Y, pch=19, col=lab, main="10% connected")
#' plot(out3$Y, pch=19, col=lab, main="25% connected")
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
do.lapeig <- function(X, ndim=2, ...){
  #------------------------------------------------------------------------
  # PREPROCESSING
  # explicit
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)

  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.lapeig : 'ndim' is a positive integer in [1,#(covariates)].")
  }

  # implicit
  params = list(...)
  pnames = names(params)

  # parameters for aux.graphnbd
  if ("type"%in%pnames){
    nbdtype = params$type
  } else {
    nbdtype = c("proportion",0.1)
  }
  if ("symmetric"%in%pnames){
    nbdsymmetric = tolower(as.character(params$symmetric))
    nbdsymmetric = match.arg(nbdsymmetric, c("union","intersect","asymmetric"))
  } else {
    nbdsymmetric = "union"
  }

  # parameters for eigenmaps
  if ("preprocess"%in%pnames){
    algpreprocess = tolower(as.character(params$preprocess))
    algpreprocess = match.arg(algpreprocess, c("null","center","scale","cscale","whiten","decorrelate"))
  } else {
    algpreprocess = "null"
  }
  if ("weighted"%in%pnames){
    wflag = as.logical(params$weighted)
  } else {
    wflag = FALSE
  }
  if ("kernelscale"%in%pnames){
    t = as.double(params$kernelscale)
    if (!is.numeric(t)||is.na(t)||(t<=0)){
      stop("* do.lapeig : 'kernelscale' is a positive real value.")
    }
    if (is.infinite(t)){
      wflag = FALSE
    }
  } else {
    t = 1.0
  }

  #------------------------------------------------------------------------
  # COMPUTE
  # data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  n = nrow(pX)
  p = ncol(pX)

  # neighborhood selection
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)

  # compute the weight matrix
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

  # compute the embedding : sorted from the bottom
  embedding = method_eigenmaps_R(W, ndim)

  # output
  result = list()
  result$Y = embedding$eigvecs[,2:(ndim+1)]
  result$eigvals = embedding$eigvals[2:(ndim+1)]
  trfinfo$algtype = "nonlinear"
  result$algorithm = "nonlinear:LAPEIG"
  result$trfinfo  = trfinfo
  return(result)
}


# auxiliary : lapeig ------------------------------------------------------
#' @keywords internal
#' @noRd
method_eigenmaps_R <- function(W, ndim){
  N = base::nrow(W)
  matToBeDec <- base::diag(N) - W/base::rowSums(W)
  eigToBeDec <- RSpectra::eigs(matToBeDec, (ndim+1), which="SM")

  output = list()
  output$eigvals = Re(eigToBeDec$values)[(ndim+1):1]
  output$eigvecs = Re(eigToBeDec$vectors)[,(ndim+1):1]
  return(output)
}
