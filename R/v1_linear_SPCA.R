#' Sparse Principal Component Analysis
#'
#' Sparse PCA (\code{do.spca}) is a variant of PCA in that each loading - or, principal
#' component - should be sparse. Instead of using generic optimization package,
#' we opt for formulating a problem as semidefinite relaxation and utilizing ADMM.
#'
#' @param X an \eqn{(n\times p)} matrix whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param mu an augmented Lagrangian parameter.
#' @param rho a regularization parameter for sparsity.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations (default: 100).}
#' \item{abstol}{absolute tolerance stopping criterion (default: 1e-8).}
#' \item{reltol}{relative tolerance stopping criterion (default: 1e-4).}
#' \item{tolerance}{stopping criterion in a Frobenius norm .}
#' }
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are principal components.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' lab   = as.factor(iris[subid,5])
#'
#' ## try different regularization parameters for sparsity
#' out1 <- do.spca(X,ndim=2,rho=0.01)
#' out2 <- do.spca(X,ndim=2,rho=1)
#' out3 <- do.spca(X,ndim=2,rho=100)
#'
#' ## embeddings for each procedure
#' Y1 <- out1$Y; Y2 <- out2$Y; Y3 <- out3$Y
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(Y1, col=lab, pch=19, main="SPCA::rho=0.01")
#' plot(Y2, col=lab, pch=19, main="SPCA::rho=1")
#' plot(Y3, col=lab, pch=19, main="SPCA::rho=100")
#' par(opar)
#'
#' @references
#' \insertRef{zou_sparse_2006}{Rdimtools}
#'
#' \insertRef{daspremont_direct_2007}{Rdimtools}
#'
#' \insertRef{ma_alternating_2013}{Rdimtools}
#'
#' @seealso \code{\link{do.pca}}
#' @author Kisung You
#' @rdname linear_SPCA
#' @concept linear_methods
#' @export
do.spca <- function(X, ndim=2, mu=1.0, rho=1.0, ...){
  #------------------------------------------------------------------------
  # Preprocessing
  if (!is.matrix(X)){stop("* do.spca : 'X' should be a matrix.")}
  myprep = ifelse(missing(preprocess), "center", match.arg(preprocess))
  myndim = min(max(1, round(ndim)), ncol(X)-1)
  mymu   = as.double(mu)
  myrho  = as.double(rho)

  # Extra parameters
  params  = list(...)
  pnames  = names(params)

  if ("abstol"%in%pnames){
    myabstol = max(.Machine$double.eps, as.double(params$abstol))
  } else {
    myabstol = 10^{-4}
  }
  if ("reltol"%in%pnames){
    myreltol = max(.Machine$double.eps, as.double(params$reltol))
  } else {
    myreltol = 10^{-2}
  }
  if ("maxiter"%in%pnames){
    myiter = max(5, round(params$maxiter))
  } else {
    myiter = 100
  }

  #------------------------------------------------------------------------
  # Version 2 update
  output = dt_spca(X, myndim, mymu, myrho, myabstol, myreltol, myiter)
  return(structure(output, class="Rdimtools"))

#
#   #------------------------------------------------------------------------
#   ## PREPROCESSING
#   #   1. data matrix
#   aux.typecheck(X)
#   n = nrow(X)
#   p = ncol(X)
#   #   2. ndim
#   ndim = as.integer(ndim)
#   if (!check_ndim(ndim,p)){
#     stop("* do.spca : 'ndim' is a positive integer in [1,#(covariates)].")
#   }
#   #   3. preprocess
#   if (missing(preprocess)){
#     algpreprocess = "center"
#   } else {
#     algpreprocess = match.arg(preprocess)
#   }
#   #   4. parameters to be passed
#   muval      = as.double(mu)
#   rhoval     = as.double(rho)
#   abstolval  = as.double(abstol)
#   reltolval  = as.double(reltol)
#   maxiterval = as.integer(maxiter)
#
#
#   #------------------------------------------------------------------------
#   ## COMPUTATION : PRELIMINARY
#   tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
#   trfinfo = tmplist$info
#   pX      = tmplist$pX
#
#   #------------------------------------------------------------------------
#   ## COMPUTATION : MAIN ITERATIVE STEP
#   #   1. compute sample covariance matrix
#   covX = stats::cov(pX)
#   #   2. compute projection and history
#   admmrun = ADMM ::admm.spca(covX, ndim, mu=muval, rho=rhoval, abstol=abstolval,
#                             reltol=reltolval, maxiter=maxiterval)
#   projection = aux.adjprojection(admmrun$basis)
#
#   #------------------------------------------------------------------------
#   ## RETURN
#   result = list()
#   result$Y = pX%*%projection
#   result$trfinfo = trfinfo
#   result$projection = projection
#   return(result)
}
