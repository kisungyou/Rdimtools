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
#' @param dtype type of potential distance transformation; \code{"log"} or \code{"sqrt"} (default: \code{"sqrt"}).
#' @param smacof a logical; \code{TRUE} to use SMACOF for Metric MDS or \code{FALSE} to use Classical MDS (default: \code{TRUE}).
#'
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations (default: 100).}
#' \item{abstol}{absolute stopping criterion for metric MDS iterations (default: 1e-8).}
#' }
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' \donttest{
#' ## load iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' lab   = as.factor(iris[,5])
#'
#' ## compare different neighborhood sizes.
#' pca2d <- do.pca(X, ndim=2)
#' phk01 <- do.phate(X, ndim=2, k=2)
#' phk02 <- do.phate(X, ndim=2, k=5)
#' phk03 <- do.phate(X, ndim=2, k=7)
#'
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2))
#' plot(pca2d$Y, col=lab, pch=19, main="PCA")
#' plot(phk01$Y, col=lab, pch=19, main="PHATE:k=2")
#' plot(phk02$Y, col=lab, pch=19, main="PHATE:k=5")
#' plot(phk03$Y, col=lab, pch=19, main="PHATE:k=7")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{moon_visualizing_2019}{Rdimtools}
#'
#' @rdname nonlinear_PHATE
#' @concept nonlinear_methods
#' @export
do.phate <- function(X, ndim=2, k=5, alpha=10, dtype=c("sqrt","log"), smacof=TRUE, ...){
  #------------------------------------------------------------------------
  # Preprocessing
  if (!is.matrix(X)){stop("* do.phate : 'X' should be a matrix.")}
  myndim = round(ndim)
  myk     = max(2, round(k))
  myalpha = max(1, as.double(alpha))
  if (missing(dtype)){
    mydtype = "sqrt"
  } else {
    mydtype = match.arg(tolower(dtype), c("sqrt","log"))
  }
  smacof  = as.logical(smacof)

  # Extra paramters
  params  = list(...)
  pnames  = names(params)

  if ("abstol"%in%pnames){
    myeps = max(.Machine$double.eps, as.double(params$abstol))
  } else {
    myeps = 1e-8
  }
  if ("maxiter"%in%pnames){
    myiter = max(5, round(params$maxiter))
  } else {
    myiter = 100
  }

  #------------------------------------------------------------------------
  # Run PHATE from 'maotai'
  distX <- stats::dist(X)
  fun_phate <- utils::getFromNamespace("hidden_PHATE","maotai")
  run_phate <- fun_phate(distX, nbdk=myk, alpha=myalpha)

  #------------------------------------------------------------------------
  # Run with RCPP for the rest
  output = dt_phate_partial(run_phate$P, myndim, mydtype, myiter, myeps, smacof)
  return(structure(output, class="Rdimtools"))
}
