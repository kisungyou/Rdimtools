#' Principal Component Analysis
#'
#' \code{do.pca} performs a classical principal component analysis (PCA) using
#' \code{RcppArmadillo} package for faster and efficient computation.
#'
#' @param X an \eqn{(n\times p)} matrix whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param cor mode of eigendecomposition. \code{FALSE} for decomposing covariance matrix,
#' and \code{TRUE} for correlation matrix.
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{vars}{a vector containing variances of projected data onto principal components.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are principal components.}
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
#' ## try covariance & correlation decomposition
#' out1 <- do.pca(X, ndim=2, cor=FALSE)
#' out2 <- do.pca(X, ndim=2, cor=TRUE)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(out1$Y, col=lab, pch=19, main="correlation decomposition")
#' plot(out2$Y, col=lab, pch=19, main="covariance decomposition")
#' par(opar)
#' }
#'
#' @author Kisung You
#' @references
#' \insertRef{pearson_liii_1901}{Rdimtools}
#'
#' @rdname linear_PCA
#' @concept linear_methods
#' @export
do.pca <- function(X, ndim=2, cor=FALSE){
  #------------------------------------------------------------------------
  # Preprocessing
  if (!is.matrix(X)){
    stop("* do.pca : 'X' should be a matrix.")
  }
  myndim = min(max(1, round(ndim)), ncol(X)-1)
  mycor  = as.logical(cor)

  #------------------------------------------------------------------------
  # Version 2 update
  output = dt_pca(X, myndim, mycor)
  output$vars = as.vector(output$vars)
  return(structure(output, class="Rdimtools"))
}



# call for later use ------------------------------------------------------
#' @keywords internal
pydo_pca <- function(myX, mydim, mycor){
  return(do.pca(myX, ndim=mydim, cor=mycor))
}
