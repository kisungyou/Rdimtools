#' Metric Multidimensional Scaling
#'
#' Metric MDS is a nonlinear method that is solved iteratively. We adopt a
#' well-known SMACOF algorithm for updates with uniform weights over all pairwise
#' distances after initializing the low-dimensional configuration via classical MDS.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations for metric MDS updates (default: 100).}
#' \item{abstol}{stopping criterion for metric MDS iterations (default: 1e-8).}
#' }
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#'
#' @examples
#' \donttest{
#' ## load iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' lab   = as.factor(iris[,5])
#'
#' ## compare with other methods
#' pca2d <- do.pca(X, ndim=2)
#' cmd2d <- do.mds(X, ndim=2)
#' mmd2d <- do.mmds(X, ndim=2)
#'
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(pca2d$Y, col=lab, pch=19, main="PCA")
#' plot(cmd2d$Y, col=lab, pch=19, main="Classical MDS")
#' plot(mmd2d$Y, col=lab, pch=19, main="Metric MDS")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{leeuw_applications_1977}{Rdimtools}
#'
#' \insertRef{borg_modern_2010}{Rdimtools}
#'
#' @rdname nonlinear_MMDS
#' @concept nonlinear_methods
#' @export
do.mmds <- function(X, ndim=2, ...){
  #------------------------------------------------------------------------
  # Preprocessing
  if (!is.matrix(X)){stop("* do.mmds : 'X' should be a matrix.")}
  p = base::ncol(X)
  myndim = as.integer(ndim)
  if (!check_ndim(myndim,p)){
    stop("* do.mmds : 'ndim' is a positive integer in [1,#(covariates)].")
  }

  # Extra parameters
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
  # Version 2 update
  output = dt_mmds(X, myndim, myiter, myeps)
  return(structure(output, class="Rdimtools"))
}
