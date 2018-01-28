#' Linear Local Tangent Space Alignment
#'
#' Linear Local Tangent Space Alignment (LLTSA) is a linear variant of the
#' celebrated LTSA method. It uses the tangent space in the neighborhood for each data point
#' to represent the local geometry. Alignment of those local tangent spaces in the low-dimensional space
#' returns an explicit mapping from the high-dimensional space.
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
#' Default is "center" and  other methods "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
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
#' ## generate data
#' X <- aux.gensamples()
#'
#' ## try different neighborhood size
#' out1 <- do.lltsa(X, type=c("proportion",0.01))
#' out2 <- do.lltsa(X, type=c("proportion",0.05))
#' out3 <- do.lltsa(X, type=c("proportion",0.10))
#'
#' ## Visualize three different projections
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1],out1$Y[,2],main="LLTSA::1% connected")
#' plot(out2$Y[,1],out2$Y[,2],main="LLTSA::5% connected")
#' plot(out3$Y[,1],out3$Y[,2],main="LLTSA::10% connected")
#' }
#'
#' @references
#' \insertRef{zhang_linear_2007}{Rdimtools}
#'
#' @seealso \code{\link{do.ltsa}}
#' @author Kisung You
#' @rdname linear_LLTSA
#' @export
do.lltsa <- function(X, ndim=2, type=c("proportion",0.1),
                     symmetric=c("union","intersect","asymmetric"),
                     preprocess=c("center","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lltsa : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. nbd setup
  nbdtype = type
  if (missing(symmetric)){
    nbdsymmetric = "asymmetric"
  } else {
    nbdsymmetric = match.arg(symmetric)
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocess of data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
  #   2. build neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION : LINEAR LTSA IN A DESCRIPTIVE ORDER
  #   1. main outer iteration
  B = array(0,c(n,n))
  for (i in 1:n){
    #   1-1. find the index
    Ii = which(nbdmask[i,])
    ki = length(Ii)
    if (ki<=1){
      stop("* do.lltsa : please choose larger neighborhood.")
    }
    #   1-2. select Xi and create Hk (centering)
    Xi = t(pX[Ii,])
    Hk = (diag(ki)-(array(1,c(ki,ki))/ki))
    #   1-3. extracting local information
    Vi = RSpectra::svds(Xi%*%Hk, ndim)$v
    Wi = Hk%*%(diag(ki)-(Vi%*%t(Vi)))
    #   1-4. alignment matrix
    B[Ii,Ii] = B[Ii,Ii] + Wi
  }
  #   2. build cost function
  LHS = (t(pX)%*%B%*%pX)
  RHS = (t(pX)%*%pX)
  #   3. extract projection matrix
  projection = aux.geigen(LHS, RHS, ndim, maximal=FALSE)

  #------------------------------------------------------------------------
  ## RETURN THE RESULTS
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
