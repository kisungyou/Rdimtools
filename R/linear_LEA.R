#' Locally Linear Embedded Eigenspace Analysis
#'
#' Locally Linear Embedding (LLE) is a powerful nonlinear manifold learning method. This method,
#' Locally Linear Embedded Eigenspace Analysis - LEA, in short - is a linear approximation to LLE,
#' similar to Neighborhood Preserving Embedding. In our implementation, the choice of weight binarization
#' is removed in order to respect original work. For 1-dimensional projection, which is rarely performed,
#' authors provided a detour for rank correcting mechanism but it is omitted for practical reason.
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
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate samples
#' X <- aux.gensamples(n=123)
#'
#' ## compare LEA with LLE and another approximation NPE
#' out1 <- do.lle(X, ndim=2)
#' out2 <- do.npe(X, ndim=2)
#' out3 <- do.lea(X, ndim=2)
#'
#' ## visual comparison
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="LLE")
#' plot(out2$Y[,1], out2$Y[,2], main="NPE")
#' plot(out3$Y[,1], out3$Y[,2], main="LEA")
#' }
#'
#' @references
#' \insertRef{fu_locally_2005}{Rdimtools}
#'
#' @seealso \code{\link{do.npe}}
#' @author Kisung You
#' @rdname linear_LEA
#' @export
do.lea <- function(X, ndim=2, type=c("proportion",0.1), symmetric=c("union","intersect","asymmetric"),
                   preprocess = c("center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lea : 'ndim' is a positive integer in [1,#(covariates)].")  }
  if (ndim==1){
    warning(" do.lea : for 1-dimensional projection, we have not implemented it yet.")
  }
  #   3. symmetric
  if (missing(symmetric)){nbdsymmetric="union"} else {nbdsymmetric=match.arg(symmetric)}
  #   4. preprocess
  if (missing(preprocess)){algpreprocess="center"} else {algpreprocess=match.arg(preprocess)}
  #   5. nbdtype
  if (missing(type)){nbdtype = c("proportion",0.1)} else {nbdtype = type}

  #------------------------------------------------------------------------
  ## COMPUTATION Part 1. Preliminary Computations
  #   1. data preprecessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. neighborhood type
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION Part 2. Main Computation for LEA
  # 1. get ready for W matrix
  W = array(0,c(n,n))
  # 2. compute rowwise elements of W
  for (i in 1:n){
    #   2-1. neighborhood index and the length
    nbdidx = which(nbdmask[i,])
    nbdnum = length(nbdidx)
    if (nbdnum<=1){
      stop("* do.lea : select larger neighborhood size. it is too less.")
    }
    #   2-2. construct Gi, solve for wi, and assign to W matrix
    W[,nbdidx] = lea_constructG_and_w(pX[i,],pX[nbdidx,])
  }
  # 3. learn embedding : smallest ones from {2 to (ndim+1)}
  IW = diag(n)-W
  LHS = (t(pX)%*%t(IW)%*%IW%*%pX)
  RHS = (t(pX)%*%pX)

  projtmp    = aux.geigen(LHS, RHS, (ndim+1), maximal=FALSE)
  projection = aux.adjprojection(projtmp[,2:(ndim+1)])


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}




#' @keywords internal
#' @noRd
lea_constructG_and_w <- function(tgtvec,tgtmat){
  k = nrow(tgtmat)
  G = array(0,c(k,k))
  for (i in 1:(k-1)){
    vec1 = as.vector(tgtvec)-as.vector(tgtmat[i,])
    for (j in (i+1):k){
      vec2 = as.vector(tgtvec)-as.vector(tgtmat[j,])
      val12 = sum(vec1*vec2)
      G[i,j] = val12
      G[j,i] = val12
    }
  }
  onesk = rep(1,k)
  term1 = solve(G,onesk)
  term2 = sum(as.vector(onesk)*as.vector(term1))
  output = term1/term2
  return(output)
}
