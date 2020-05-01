#' Curvilinear Distance Analysis
#'
#' Curvilinear Distance Analysis (CRDA) is a variant of Curvilinear Component Analysis in that
#' the input pairwise distance is altered by curvilinear distance on a data manifold.
#' Like in Isomap, it first generates \emph{neighborhood graph} and finds \emph{shortest path} on
#' a constructed graph so that the shortest-path length plays as an approximate geodesic distance on
#' nonlinear manifolds.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#' \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#' Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#' among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.
#' @param weight \code{TRUE} to perform CRDA on weighted graph, or \code{FALSE} otherwise.
#' @param lambda threshold value.
#' @param alpha initial value for updating.
#' @param maxiter maximum number of iterations allowed.
#' @param tolerance stopping criterion for maximum absolute discrepancy between two distance matrices.
#'
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{niter}{the number of iterations until convergence.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#'
#' @examples
#' ## load iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' ## different settings of connectivity
#' out1 <- do.crda(X, type=c("proportion",0.10))
#' out2 <- do.crda(X, type=c("proportion",0.25))
#' out3 <- do.crda(X, type=c("proportion",0.50))
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, pch=19, main="10% connected")
#' plot(out2$Y, col=label, pch=19, main="25% connected")
#' plot(out3$Y, col=label, pch=19, main="50% connected")
#' par(opar)
#'
#' @references
#' \insertRef{lee_curvilinear_2002}{Rdimtools}
#'
#' \insertRef{lee_nonlinear_2004}{Rdimtools}
#'
#' @seealso \code{\link{do.isomap}}, \code{\link{do.crca}}
#' @author Kisung You
#' @rdname nonlinear_CRDA
#' @concept nonlinear_methods
#' @export
do.crda <- function(X,ndim=2,type=c("proportion",0.1),symmetric="union",weight=TRUE,
                    lambda=1.0,alpha=1.0,maxiter=1000,tolerance=1e-6){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  vlnum = 1e+10
  # 1. data X
  aux.typecheck(X)
  # 2. ndim
  ndim = as.integer(ndim)
  if(!check_NumMM(ndim,1,ncol(X),compact=FALSE)){stop("* do.crda : 'ndim' is a positive integer in (1,#(covariates)).")}
  # 3. lambda
  lambda = as.double(lambda)
  if(!check_NumMM(lambda,0,vlnum,compact=FALSE)){stop("* do.crda : 'lambda' should be a positive real number.")}
  # 4. alpha
  alpha = as.double(alpha)
  if(!check_NumMM(alpha,0,vlnum,compact=FALSE)){stop("* do.crda : 'alpha' should be a positive real number.")}
  # 5. maxiter
  maxiter = as.integer(maxiter)
  if(!check_NumMM(maxiter,1,vlnum,compact=FALSE)){stop("* do.crda : 'maxiter' should be a not-so-small positive integer.")}
  # 6. tolerance
  tolerance = as.double(tolerance)
  if(!check_NumMM(alpha,0,vlnum,compact=FALSE)){stop("* do.crda : 'tolerance' should be a positive real number.")}
  # 7. ISOMAP parameters
  nbdtype = type
  nbdsymmetric = symmetric
  if (!is.element(nbdsymmetric,c("union","intersect","asymmetric"))){
    stop("* do.crda : 'symmetric' should have one of three values.")
  }
  algweight = weight
  if (!is.logical(algweight)){
    stop("* do.crda : 'weight' should be a logical value.")
  }

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. ISOMAP-type Curvilinear Distance
  nbdstruct = aux.graphnbd(X,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  D     = nbdstruct$dist
  Dmask = nbdstruct$mask
  nD    = ncol(D)
  if (algweight){
    wD = Dmask*D
    idnan = is.na(wD)
    wD[idnan] = 0
  } else {
    wD = matrix(as.double(Dmask),nrow=nD)
  }
  Xij = aux.shortestpath(wD)
  #   2. Initialization via PCA
  Yinit = do.pca(X, ndim=ndim)$Y
  #   3. vecselector for random-number generation
  vecselector = as.vector(sample(0:(nrow(X)-1), maxiter, replace=TRUE))
  #   4. main computation
  Youtput = method_crca(Xij,Yinit,lambda,alpha,maxiter,tolerance,vecselector)

  #------------------------------------------------------------------------
  ## RETURN OUTPUT
  trfinfo = list()
  trfinfo$type = "null"
  trfinfo$algtype = "nonlinear"
  result = list()
  result$Y = Youtput$Y
  result$niter = Youtput$niter
  result$trfinfo = trfinfo
  return(result)
}


