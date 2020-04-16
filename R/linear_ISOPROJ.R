#' Isometric Projection
#'
#' Isometric Projection is a linear dimensionality reduction algorithm that exploits
#' geodesic distance in original data dimension and mimicks the behavior in the target dimension.
#' Embedded manifold is approximated by graph construction as of ISOMAP. Since it involves
#' singular value decomposition and guesses intrinsic dimension by the number of positive singular values
#' from the decomposition of data matrix, it automatically corrects the target dimension accordingly.
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
#' \item{Y}{an \eqn{(n\times ndim)} matrix of projected observations as rows.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are loadings.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \donttest{
#' ## generate data
#' X <- aux.gensamples(n=28)
#'
#' ## 1. connecting 10% of data for graph construction.
#' output1 <- do.isoproj(X,ndim=2,type=c("proportion",0.10))
#'
#' ## 2. constructing 25%-connected graph
#' output2 <- do.isoproj(X,ndim=2,type=c("proportion",0.25))
#'
#' ## Visualize two different projections
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(output1$Y, main="10%")
#' plot(output2$Y, main="25%")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{cai_isometric_2007}{Rdimtools}
#'
#' @rdname linear_ISOPROJ
#' @author Kisung You
#' @concept linear_methods 
#' @export
do.isoproj <- function(X,ndim=2,type=c("proportion",0.1),
                       symmetric=c("union","intersect","asymmetric"),
                       preprocess=c("center","scale","cscale","decorrelate","whiten")){
  ## PREPROCESSING
  # Preprocessing : typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("*do.isoproj : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # Parameters
  # 1. aux.graphnbd
  #   type : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #   symmetric : 'intersect','union', or 'asymmetric'
  # 2. isomap itself
  #   preprocess : 'center','decorrelate', or 'whiten'
  nbdtype = type
  if (missing(symmetric)){
    nbdsymmetric = "union"
  } else {
    nbdsymmetric = match.arg(symmetric)
  }

  algweight = TRUE
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  ## COMPUTATION
  #   1. data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. neighborhood selection
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  D     = nbdstruct$dist
  Dmask = nbdstruct$mask
  nD    = ncol(D)
  #   3. network binarization
  if (algweight){
    wD = Dmask*D
    idnan = is.na(wD)
    wD[idnan] = 0
  } else {
    wD = matrix(as.double(Dmask),nrow=nD)
  }
  #   4. compute shortest path
  sD = aux.shortestpath(wD)
  #   5. convert into gram matrix
  rDG = convert_gram(sD)
  #   6. SVD using RSpectra : if rank-deficient, change!
  workdim = RSpectra::svds(t(pX), ndim)
  if (workdim$d[ndim] <= 0){
    message("do.isoproj : effective rank of input matrix X is less than 'ndim'. Change 'ndim' correspondingly.")
    ndim = sum((workdim$d>0))
  }
  U = workdim$u[,1:ndim]
  V = workdim$v[,1:ndim]
  Xbar = t(pX%*%U)  # now it corresponds to original notation.
  #   7. generalized eigenvalue problem
  LHS = Xbar%*%rDG%*%t(Xbar)
  RHS = Xbar%*%t(Xbar)
  SOL = aux.bicgstab(RHS, LHS, verbose=FALSE)
  #   8. eigendecomposition and compute projection matrix A
  if (ndim<3){
    res = base::eigen(SOL$x, ndim)
    A   = U%*%(res$vectors[,1:ndim])
  } else {
    res = RSpectra::eigs_sym(SOL$x, ndim, which="LA")
    A = U%*%(res$vectors)
  }


  ## RETURN OUTPUT
  projection = aux.adjprojection(A)
  result = list()
  result$Y = pX%*%projection
  result$projection = projection
  result$trfinfo    = trfinfo
  return(result)
}


#' @keywords internal
#' @noRd
convert_gram <- function(D){
  m = nrow(D)
  H = diag(x=1,m,m)-(outer(rep(1,m),rep(1,m))/m)
  S = (D^2)
  return((-H%*%S%*%H)/2)
}
