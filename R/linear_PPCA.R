#' Probabilistic Principal Component Analysis
#'
#' Probabilistic PCA (PPCA) is a probabilistic framework to explain the well-known PCA model. Using
#' the conjugacy of normal model, we compute MLE for values explicitly derived in the paper. Note that
#' unlike PCA where loadings are directly used for projection, PPCA uses \eqn{WM^{-1}} as projection matrix,
#' as it is relevant to the error model. Also, for high-dimensional problem, it is possible that MLE can have
#' negative values if sample covariance given the data is rank-deficient.
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{mle.sigma2}{MLE for \eqn{\sigma^2}.}
#' \item{mle.W}{MLE of a \eqn{(p\times ndim)} mapping from latent to observation in column major.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150, 50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' ## Compare PCA and PPCA
#' PCA  <- do.pca(X, ndim=2)
#' PPCA <- do.ppca(X, ndim=2)
#'
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(PCA$Y,  pch=19, col=label, main="PCA")
#' plot(PPCA$Y, pch=19, col=label, main="PPCA")
#' par(opar)
#' }
#'
#' @seealso \code{\link{do.pca}}
#' @author Kisung You
#' @references
#' \insertRef{tipping_probabilistic_1999}{Rdimtools}
#'
#' @rdname linear_PPCA
#' @concept linear_methods
#' @export
do.ppca <- function(X, ndim=2){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  # 1. data X
  aux.typecheck(X)
  # 2. ndim
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>=ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("*do.pca : 'ndim' is a positive integer in [1,#(covariates)).")
  }
  # # 3. preprocess
  # if (missing(preprocess)){
  #   algpreprocess = "center"
  # } else {
  #   algpreprocess = match.arg(preprocess)
  # }

  #------------------------------------------------------------------------
  ## COMPUTATION
  # #   1. Preprocessing the data
  # tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  # trfinfo = tmplist$info
  # matT    = tmplist$pX
  #   2. parameter
  d = ncol(X)
  q = as.integer(ndim)
  #   3. sample covariance and eigendecomposition
  S    = stats::cov(X)
  eigS = base::eigen(S, TRUE)
  #   4. ML estimates
  mlsig2 = (sum(eigS$values[(q+1):d]))/(d-q)
  mlW    = (eigS$vectors[,1:q])%*%(diag(eigS$values[1:q] - mlsig2))
  #   5. Projection
  M = (t(mlW)%*%mlW)+(diag(ncol(mlW))*mlsig2)
  SOL = aux.bicgstab(M, t(mlW), verbose=FALSE)
  projection = aux.adjprojection(t(SOL$x))

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y          = (X%*%projection)
  result$projection = projection
  result$mle.sigma2  = mlsig2
  result$mle.W       = mlW
  result$algorithm   = "linear:PPCA"
  return(structure(result, class="Rdimtools"))
}





