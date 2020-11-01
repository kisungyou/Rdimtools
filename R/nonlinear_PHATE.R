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
#' plot(phlog$Y, col=lab, pch=19, main="'log' distance")
#' plot(phsqt$Y, col=lab, pch=19, main="'sqrt' distance")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{moon_visualizing_2019}{Rdimtools}
#'
#' @rdname nonlinear_PHATE
#' @concept nonlinear_methods
#' @export
do.phate <- function(X, ndim=2, k=5, alpha=10, dtype=c("log","sqrt"), maxiter=100, abstol=1e-6){
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


  #------------------------------------------------------------------------
  # Version 2 update
  output = dt_phate(X, myndim, myprep, myk, myalpha, mydtype, myiter, myeps)
  return(output)
}

# vecnorm <- function(x){
#   return(sqrt(sum(x^2)))
# }

# X = as.matrix(iris[sample(1:150, 30),1:4])
# D = as.matrix(dist(X))
# K    = exp(-D^2)
# Ksum = rowSums(K)
#
# P = diag(1/Ksum)%*%K
# A = diag(1/sqrt(Ksum))%*%K%*%diag(1/sqrt(Ksum))
#
# eigA = eigen(A)$values
#
# A2 = A%*%A
# A3 = A2%*%A
# A4 = A3%*%A
# A5 = A4%*%A
#
# vecnorm(eigen(A2)$values - eigA^2)
# vecnorm(eigen(A3)$values - eigA^3)
# vecnorm(eigen(A4)$values - eigA^4)
# vecnorm(eigen(A5)$values - eigA^5)

