#' Locality-Preserved Maximum Information Projection
#'
#' Locality-Preserved Maximum Information Projection (LPMIP) is an unsupervised linear dimension reduction method
#' to identify the underlying manifold structure by learning both the within- and between-locality information. The
#' parameter \code{alpha} is balancing the tradeoff between two and the flexibility of this model enables an interpretation
#' of it as a generalized extension of LPP.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null" and three options of "center", "decorrelate", or "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param sigma bandwidth parameter for heat kernel in \eqn{(0,\infty)}.
#' @param alpha balancing parameter between two locality information in \eqn{[0,1]}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## generate default dataset
#' X <- aux.gensamples(n=123)
#'
#' ## try different neighborhood size
#' out1 <- do.lpmip(X, ndim=2, type=c("proportion",0.01))
#' out2 <- do.lpmip(X, ndim=2, type=c("proportion",0.1))
#' out3 <- do.lpmip(X, ndim=2, type=c("proportion",0.25))
#'
#' ## Visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1],out1$Y[,2],main="1% connected")
#' plot(out2$Y[,1],out2$Y[,2],main="10% connected")
#' plot(out3$Y[,1],out3$Y[,2],main="25% connected")
#'
#' @references
#' \insertRef{haixian_wang_locality-preserved_2008}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LPMIP
#' @export
do.lpmip <- function(X, ndim=2, type=c("proportion",0.1), preprocess=c("center","whiten","decorrelate"),
                     sigma=10, alpha=0.5){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lpmip : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. neighborhood information
  nbdtype = type
  nbdsymmetric = "union"
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. sigma
  sigma = as.double(sigma)
  if (!check_NumMM(sigma, 0, Inf, compact=FALSE)){stop("* do.lpmip : 'sigma' is a bandwidth parameter in (0,Inf).")}
  #   6. alpha
  alpha = as.double(alpha)
  if (!check_NumMM(alpha,0,1,compact=TRUE)){stop(" do.lpmip : 'alpha' is a balancing parameter in [0,1].")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"

  #   2. build neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR LPMIP
  #   1. W : weight matrix & Ltilde
  Dsqmat = (as.matrix(dist(pX))^2)
  W      = exp(-Dsqmat/sigma)
  Ltilde = diag(rowSums(W))-W
  #   2. A with neighborhood
  A      = W*nbdmask; diag(A)=0;
  L      = diag(rowSums(A))-A
  #   3. cost function
  costW = t(pX)%*%(alpha*Ltilde - L)%*%pX
  #   4. compute projection vectors
  projection = aux.adjprojection(RSpectra::eigs(costW, ndim)$vectors)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
