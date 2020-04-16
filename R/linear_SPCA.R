#' Sparse Principal Component Analysis
#'
#' Sparse PCA (\code{do.spca}) is a variant of PCA in that each loading - or, principal
#' component - should be sparse. Instead of using generic optimization package,
#' we opt for formulating a problem as semidefinite relaxation and utilizing ADMM.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param mu an augmented Lagrangian parameter.
#' @param rho a regularization parameter for sparsity.
#' @param abstol absolute tolerance stopping criterion.
#' @param reltol relative tolerance stopping criterion.
#' @param maxiter maximum number of iterations.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are principal components.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{history}{a length-\code{ndim} list where each element is an iteration history.}
#' }
#'
#'
#' @examples
#' \donttest{
#' ## generate default dataset and make its dimension three-folds.
#' Xp <- aux.gensamples()
#' X  <- cbind(Xp,Xp,Xp)
#'
#' ## try different regularization parameters for sparsity
#' out1 <- do.spca(X,ndim=2,rho=0.01)
#' out2 <- do.spca(X,ndim=2,rho=1)
#' out3 <- do.spca(X,ndim=2,rho=100)
#'
#' ## Visualize principal components as columns in an image
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' image(t(out1$projection), main="SPCA::rho=0.01")
#' image(t(out2$projection), main="SPCA::rho=1")
#' image(t(out3$projection), main="SPCA::rho=100")
#' par(opar)
#' }
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
do.spca <- function(X, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"),
                    mu=1.0, rho=1.0, abstol=1e-4, reltol=1e-2, maxiter=1000){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.spca : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. parameters to be passed
  muval      = as.double(mu)
  rhoval     = as.double(rho)
  abstolval  = as.double(abstol)
  reltolval  = as.double(reltol)
  maxiterval = as.integer(maxiter)


  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN ITERATIVE STEP
  #   1. compute sample covariance matrix
  covX = cov(pX)
  #   2. compute projection and history
  admmrun = ADMM::admm.spca(covX, ndim, mu=muval, rho=rhoval, abstol=abstolval,
                            reltol=reltolval, maxiter=maxiterval)
  projection = aux.adjprojection(admmrun$basis)
  history    = admmrun$history

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  result$history = history
  return(result)
}
