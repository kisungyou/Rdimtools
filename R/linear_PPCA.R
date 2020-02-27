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
#' @param preprocess an option for preprocessing the data. This supports three methods.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are principal components.}
#' \item{mle.sigma2}{MLE for \eqn{\sigma^2}.}
#' \item{mle.W}{MLE of a \eqn{(p\times ndim)} mapping from latent to observation in column major.}
#' }
#'
#' @examples
#' \dontrun{
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.integer(iris$Species)
#'
#' ## Compare PCA and PPCA
#' PCA  <- do.pca(X, ndim=2, preprocess="center")
#' PPCA <- do.ppca(X, ndim=2, preprocess="center")
#'
#' ## Visualize
#' opar <- par(mfrow=c(1,2))
#' plot(PCA$Y,  col=label, main="PCA")
#' plot(PPCA$Y, col=label, main="PPCA")
#' par(opar)
#' }
#'
#' @seealso \code{\link{do.pca}}
#' @author Kisung You
#' @references
#' \insertRef{tipping_probabilistic_1999}{Rdimtools}
#'
#' @rdname linear_PPCA
#' @export
do.ppca <- function(X, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  # 1. data X
  aux.typecheck(X)
  # 2. ndim
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>=ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("*do.pca : 'ndim' is a positive integer in [1,#(covariates)).")
  }
  # 3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION
  #   1. Preprocessing the data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  matT    = tmplist$pX
  #   2. parameter
  d = ncol(matT)
  q = as.integer(ndim)
  #   3. sample covariance and eigendecomposition
  S    = cov(matT)
  eigS = base::eigen(S, TRUE)
  #   4. ML estimates
  mlsig2 = (sum(eigS$values[(q+1):d]))/(d-q)
  mlW    = (eigS$vectors[,1:q])%*%(diag(eigS$values[1:q] - mlsig2))
  #   5. Projection
  M = (t(mlW)%*%mlW)+(diag(ncol(mlW))*mlsig2)
  SOL = aux.bicgstab(M, t(mlW), verbose=FALSE)
  projection = t(SOL$x)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y          = (matT%*%projection)
  result$trfinfo    = trfinfo
  result$projection = projection
  result$mle.sigma2  = mlsig2
  result$mle.W       = mlW
  return(result)
}





