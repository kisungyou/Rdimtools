#' Local Tangent Space Alignment
#'
#' Local Tangent Space Alignment, or LTSA in short, is a nonlinear dimensionality reduction method
#' that mimicks the behavior of low-dimensional manifold embedded in high-dimensional space.
#' Similar to LLE, LTSA computes tangent space using nearest neighbors of a given data point, and
#' a multiple of tangent spaces are gathered to to find an embedding that aligns the tangent spaces
#' in target dimensional space.
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
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{eigvals}{a vector of eigenvalues from the final decomposition.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate data
#' X <- aux.gensamples(dname="cswiss",n=100)
#'
#' ## 1. use 10%-connected graph
#' output1 <- do.ltsa(X,ndim=2)
#'
#' ## 2. use 25%-connected graph
#' output2 <- do.ltsa(X,ndim=2,type=c("proportion",0.25))
#'
#' ## 3. use 50%-connected graph
#' output3 <- do.ltsa(X,ndim=2,type=c("proportion",0.50))
#'
#' ## Visualize three different projections
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(output1$Y, main="10%")
#' plot(output2$Y, main="25%")
#' plot(output3$Y, main="50%")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zhang_linear_2007}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_LTSA
#' @export
do.ltsa <- function(X, ndim=2, type=c("proportion",0.1),
                    symmetric=c("union","intersect","asymmetric"),
                    preprocess=c("center","scale","cscale","decorrelate","whiten")){
  # process : typechecking
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.ltsa : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  nbdtype = type
  if (missing(symmetric)){
    nbdsymmetric = "union"
  } else {
    nbdsymmetric = match.arg(symmetric)
  }
  algweight = FALSE
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  # process : data processing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  ## Main Computations
  # 1. neighborhood selection
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask
  # 2. iterative steps, let's go
  K = array(0,c(nrow(nbdmask),ncol(nbdmask)))
  for (i in 1:nrow(nbdmask)){
    # 2-1. local coordinate relations
    idxNi = nbdmask[1,]
    Xi    = (scale(pX[idxNi,],center=TRUE,scale=FALSE))
    ki    = nrow(Xi)
    # 2-2. svd, Vi and Gi
    outsvd = eigen(Xi%*%t(Xi))
    Vi     = outsvd$vectors[,(1:min(ndim,ncol(outsvd$v)))] # dim : (k-by-d)
    Gi     = cbind(matrix(1,c(ki,1))/ki, Vi)
    Wi     = diag(ki) - (Gi%*%t(Gi))
    # 2-3. LTSA kernel construction via global alignment
    K[idxNi,idxNi] = K[idxNi,idxNi] + Wi
  }
  # 2-4. force B to be symmetric
  K = (K+t(K))/2

  # 3. eigendecomposition
  #    i think we need to find largest
  eigK   = eigen(K) # in a decreasing order for 'values'
  seqvec = order(eigK$values,decreasing=FALSE)[1:(ndim+1)]

  if (min(eigK$values)==0){
    outval = eigK$values[seqvec[2:(ndim+1)]]
    outvec = eigK$vectors[,seqvec[2:(ndim+1)]]
  } else {
    outval = eigK$values[seqvec[1:ndim]]
    outvec = eigK$vectors[,seqvec[1:ndim]]
  }

  ## Return output
  result = list()
  result$Y = outvec
  result$trfinfo = trfinfo
  result$eigvals = outval
  return(result)
}
