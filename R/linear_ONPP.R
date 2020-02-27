#' Orthogonal Neighborhood Preserving Projections
#'
#' Orthogonal Neighborhood Preserving Projection (ONPP) is an unsupervised linear dimension reduction method.
#' It constructs a weighted data graph from LLE method. Also, it develops LPP method by preserving
#' the structure of local neighborhoods.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
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
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.integer(iris$Species)
#'
#' ## try different numbers for neighborhood size
#' out1 = do.onpp(X, type=c("proportion",0.05))
#' out2 = do.onpp(X, type=c("proportion",0.1))
#' out3 = do.onpp(X, type=c("proportion",0.25))
#'
#' ## visualize
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(out1$Y, col=label, main="ONPP::5% connectivity")
#' plot(out2$Y, col=label, main="ONPP::10% connectivity")
#' plot(out3$Y, col=label, main="ONPP::25% connectivity")
#' par(opar)
#'
#' @references
#' \insertRef{kokiopoulou_orthogonal_2007}{Rdimtools}
#'
#' @rdname linear_ONPP
#' @author Kisung You
#' @export
do.onpp <- function(X, ndim=2, type=c("proportion",0.1),
                    preprocess=c("center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.onpp : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   3. type
  nbdtype = type
  nbdsymmetric = "union"
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

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
  ## COMPUTATION : MAIN COMPUTATION FOR ONPP
  #   1. cost function : Mtilde
  diagN  = diag(n)
  Mtilde = (t(pX)%*%t(diagN-W)%*%(diagN-W)%*%pX)

  #   2. use (ndim+1) lowest, except for the first one
  if ((ndim+1)==p){
    eigM = base::eigen(Mtilde)
    projection = aux.adjprojection(eigM$vectors[,p:2])
  } else {
    projection = aux.adjprojection(RSpectra::eigs(Mtilde, ndim+1, which="SR")$vectors[,2:(ndim+1)])
  }

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
