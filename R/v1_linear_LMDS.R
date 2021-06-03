#' Landmark Multidimensional Scaling
#'
#' Landmark MDS is a variant of Classical Multidimensional Scaling in that
#' it first finds a low-dimensional embedding using a small portion of given dataset
#' and graft the others in a manner to preserve as much pairwise distance from
#' all the other data points to landmark points as possible.
#'
#' @param X an \eqn{(n\times p)} matrix whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param npoints the number of landmark points to be drawn.
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' lab   = as.factor(iris[,5])
#'
#' ## use 10% and 25% of the data and compare with full MDS
#' output1 <- do.lmds(X, ndim=2, npoints=round(nrow(X)*0.10))
#' output2 <- do.lmds(X, ndim=2, npoints=round(nrow(X)*0.25))
#' output3 <- do.mds(X, ndim=2)
#'
#' ## vsualization
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(output1$Y, pch=19, col=lab, main="10% random points")
#' plot(output2$Y, pch=19, col=lab, main="25% random points")
#' plot(output3$Y, pch=19, col=lab, main="original MDS")
#' par(opar)
#' }
#'
#' @seealso \code{\link{do.mds}}
#' @references
#' \insertRef{silva_global_2002}{Rdimtools}
#'
#' \insertRef{lee_landmark_2009}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LMDS
#' @concept linear_methods
#' @export
do.lmds <- function(X, ndim=2, npoints=max(nrow(X)/5,ndim+1)){
  #------------------------------------------------------------------------
  # Preprocessing
  if (!is.matrix(X)){stop("* do.lmds : 'X' should be a matrix.")}
  myndim = min(max(1, round(ndim)), ncol(X)-1)
  mynpts = max(round(npoints),2)

  #------------------------------------------------------------------------
  # Version 2 update
  output = dt_lmds(X, myndim, mynpts)
  return(structure(output, class="Rdimtools"))
}


# call for later use ------------------------------------------------------
#' @keywords internal
pydo_lmds <- function(myX, mydim, mynpts){
  return(do.lmds(myX, ndim=mydim, npoints=mynpts))
}
