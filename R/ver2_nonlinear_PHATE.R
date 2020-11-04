#' Potential of Heat Diffusion for Affinity-based Transition Embedding
#'
#' PHATE is a nonlinear method that is specifically targeted at visualizing
#' high-dimensional data by embedding it on 2- or 3-dimensional space. We offer
#' a native implementation of PHATE solely in R/C++ without interface to python module.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param k size of nearest neighborhood (default: 5).
#' @param alpha decay parameter for Gaussian kernel exponent (default: 10).
#' @param dtype type of potential distance transformation; \code{"log"} or \code{"sqrt"}.
#' @param smacof a logical; \code{TRUE} to use SMACOF for Metric MDS (default), or \code{FALSE} to use Classical MDS.
#' @param maxiter maximum number of iterations for metric MDS updates (default: 100).
#' @param abstol stopping criterion for metric MDS iterations (default: 1e-6).
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \donttest{
#' ## load iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' lab   = as.factor(iris[,5])
#'
#' ## compare two potential distances with PCA
#' pca2d <- do.pca(X, ndim=2)
#' phlog <- do.phate(X, ndim=2, dtype="log")
#' phsqt <- do.phate(X, ndim=2, dtype="sqrt")
#'
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(pca2d$Y, col=lab, pch=19, main="PCA")
#' plot(phlog$Y, col=lab, pch=19, main="log potential")
#' plot(phsqt$Y, col=lab, pch=19, main="sqrt potential")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{moon_visualizing_2019}{Rdimtools}
#'
#' @rdname nonlinear_PHATE
#' @concept nonlinear_methods
#' @export
do.phate <- function(X, ndim=2, k=5, alpha=10, dtype=c("log","sqrt"), smacof=TRUE, maxiter=100, abstol=1e-6){
  #------------------------------------------------------------------------
  # Preprocessing
  if (!is.matrix(X)){stop("* do.phate : 'X' should be a matrix.")}
  myndim = round(ndim)
  myprep = "null"
  myk     = max(2, round(k))
  myalpha = max(1, as.double(alpha))
  mydtype = match.arg(tolower(dtype),c("log","sqrt"))
  myiter  = max(50, round(maxiter))
  myeps   = max(as.double(abstol), 100*.Machine$double.eps)
  smacof  = as.logical(smacof)


  #------------------------------------------------------------------------
  # Version 2 update - tried knn approximate search but somehow failed for symmetry
  # nbx    = nabor::knn(X, k=(myk+1))
  # nbdist = as.vector(nbx$nn.dists[,(myk+1)])
  output = dt_phate(X, myndim, myprep, myk, myalpha, mydtype, myiter, myeps, smacof)
  return(output)
}
