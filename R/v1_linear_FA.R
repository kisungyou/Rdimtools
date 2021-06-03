#' Exploratory Factor Analysis
#'
#' \code{do.fa} is an optimization-based implementation of a popular technique for Exploratory Data Analysis.
#' It is closely related to principal component analysis.
#'
#' @param X an \eqn{(n\times p)} matrix whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued number of loading variables, or target dimension.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations (default: 10).}
#' \item{tolerance}{stopping criterion in a Frobenius norm (default: 1e-8).}
#' }
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{loadings}{a \eqn{(p\times ndim)} matrix whose rows are extracted loading factors.}
#' \item{noise}{a length-\eqn{p} vector of estimated noise.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' lab   = as.factor(iris[subid,5])
#'
#' ## compare with PCA and MDS
#' out1 <- do.fa(X, ndim=2)
#' out2 <- do.mds(X, ndim=2)
#' out3 <- do.pca(X, ndim=2)
#'
#' ## visualize three different projections
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=lab, main="Factor Analysis")
#' plot(out2$Y, pch=19, col=lab, main="MDS")
#' plot(out3$Y, pch=19, col=lab, main="PCA")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{spearman_general_1904}{Rdimtools}
#'
#' @rdname linear_FA
#' @author Kisung You
#' @concept linear_methods
#' @export
do.fa <- function(X, ndim=2, ...){
  #------------------------------------------------------------------------
  # Basic
  if (!is.matrix(X)){
    stop("* do.fa : 'X' should be a matrix.")
  }
  myndim = min(max(1, round(ndim)), ncol(X)-1)

  # Extra Parameters
  params  = list(...)
  pnames  = names(params)

  if ("maxiter" %in% pnames){
    myiter = max(5, round(params$maxiter))
  } else {
    myiter = 10
  }
  if (("tolerance"%in%pnames)){
    mytol = max(.Machine$double.eps, as.double(params$tolerance))
  } else {
    mytol = sqrt(.Machine$double.eps)
  }

  #------------------------------------------------------------------------
  # Version 2 update
  output = dt_fa(X, myndim, myiter, mytol)
  output$noise = as.vector(output$noise)
  return(structure(output, class="Rdimtools"))
}

# call for later use ------------------------------------------------------
#' @keywords internal
pydo_fa <- function(myX, mydim, myiter, myeps){
  return(do.fa(myX, ndim=mydim, maxiter=myiter, tolerance=myeps))
}
