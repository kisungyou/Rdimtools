#' Metric Multidimensional Scaling
#'
#' Metric MDS is a nonlinear method that is solved iteratively. We adopt a
#' well-known SMACOF algorithm for updates with uniform weights over all pairwise
#' distances after initializing the low-dimensional configuration via classical MDS.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param maxiter maximum number of iterations for metric MDS updates (default: 100).
#' @param abstol stopping criterion for metric MDS iterations (default: 1e-6).
#' @param preprocess an option for preprocessing the data. Default is \code{"null"}.
#' See also \code{\link{aux.preprocess}} for more details.
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
do.mmds <- function(X,ndim=2,
                    preprocess=c("null","center","scale","cscale","decorrelate","whiten"),
                    maxiter=100, abstol=1e-6){
  #------------------------------------------------------------------------
  # Preprocessing
  if (!is.matrix(X)){stop("* do.mmds : 'X' should be a matrix.")}
  myndim = round(ndim)
  myprep = ifelse(missing(preprocess), "null", match.arg(preprocess))
  myiter  = max(50, round(maxiter))
  myeps   = max(as.double(abstol), 100*.Machine$double.eps)

  #------------------------------------------------------------------------
  # Version 2 update
  output = dt_mmds(X, myndim, myprep, myiter, myeps)
  return(output)
}
