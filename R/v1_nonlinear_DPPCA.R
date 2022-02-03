#' Dual Probabilistic Principal Component Analysis
#'
#' Dual view of PPCA optimizes the latent variables directly from a simple
#' Bayesian approach to model the noise using the multivariate Gaussian distribution
#' of zero mean and spherical covariance \eqn{\beta^{-1} I}. When \eqn{\beta} is too small,
#' the algorithm automatically returns an error and provides a guideline for minimal
#' value that enables successful computation.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param beta the degree for modeling the level of noise (default: 1).
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
#' ## compare difference choices of 'beta'
#' embed1 <- do.dppca(X, beta=0.2)
#' embed2 <- do.dppca(X, beta=1)
#' embed3 <- do.dppca(X, beta=5)
#'
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(embed1$Y , col=lab, pch=19, main="beta=0.2")
#' plot(embed2$Y , col=lab, pch=19, main="beta=1")
#' plot(embed3$Y , col=lab, pch=19, main="beta=5")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{lawrence_probabilistic_2005}{Rdimtools}
#'
#' @seealso \code{\link{do.ppca}}
#' @rdname nonlinear_DPPCA
#' @concept nonlinear_methods
#' @export
do.dppca <- function(X, ndim=2, beta=1){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)

  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.dppca : 'ndim' is a positive integer in [1,#(covariates)].")
  }

  #   3. 'beta' parameter
  beta = as.double(beta)
  if (beta <= .Machine$double.eps){
    stop("* do.dppca : 'beta' should be a nonnegative number.")
  }


  #------------------------------------------------------------------------
  ## COMPUTATION
  #   1. centering
  Y = as.matrix(base::scale(X, center=TRUE, scale=FALSE))

  #   2. eigen decomposition
  tgtY = (Y%*%t(Y))/p
  eigY = RSpectra::eigs_sym(tgtY, ndim, which="LM")
  print(eigY$values)

  #   3. embedding
  if (any(eigY$values < (1/beta))){
    stop(paste0("* do.dppca : try a larger beta value than ", round(max(1/eigY$values),5), "."))
  }
  L = 1/sqrt((eigY$values - 1/(beta)))

  #------------------------------------------------------------------------
  # Return
  result = list()
  result$Y = eigY$vectors%*%diag(L)
  result$algorithm = "nonlinear:DPPCA"
  return(structure(result, class="Rdimtools"))

}
