#' Principal Component Analysis
#'
#' \code{do.pca} performs a classical principal component analysis \insertCite{pearson_liii_1901}{Rdimtools} using
#' \code{RcppArmadillo} package for faster and efficient computation.
#'
#' @param X an \eqn{(n\times p)} matrix whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param ... extra parameters including \describe{
#' \item{cor}{mode of eigendecomposition. \code{FALSE} for decomposing covariance matrix (default),
#' and \code{TRUE} for correlation matrix.}
#' \item{preprocess}{an additional option for preprocessing the data.
#' Default is \code{"center"}. See also \code{\link{aux.preprocess}} for more details.}
#' }
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{vars}{a vector containing variances of projected data onto principal components.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
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
#' \insertAllCited{}
#'
#' @rdname linear_PCA
#' @concept linear_methods
#' @export
do.pca <- function(X, ndim=2, ...){
  #------------------------------------------------------------------------
  # BASIC
  # explicit
  if (!is.matrix(X)){
    stop("* do.pca : 'X' should be a matrix.")
  }
  myndim = min(max(1, round(ndim)), ncol(X)-1)

  # implicit
  params = list(...)
  pnames = names(params)

  if ("cor"%in%pnames){
    par_cor = as.logical(params$cor)
  } else {
    par_cor = FALSE
  }
  if ("preprocess"%in%pnames){
    par_preprocess = tolower(params$preprocess)
  } else {
    par_preprocess = "center"
  }

  #------------------------------------------------------------------------
  # COMPUTE
  # preprocess
  tmplist = aux.preprocess.hidden(X, type=par_preprocess, algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  # main part
  output = dt_pca(pX, myndim, par_cor)

  #------------------------------------------------------------------------
  # RETURN
  output$vars = as.vector(output$vars)
  output$trfinfo = trfinfo
  return(structure(output, class="Rdimtools"))
}



# call for later use ------------------------------------------------------
#' @keywords internal
pydo_pca <- function(myX, mydim, mycor){
  return(do.pca(myX, ndim=mydim, cor=mycor))
}
