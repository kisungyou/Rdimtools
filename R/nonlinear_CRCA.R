#' Curvilinear Component Analysis
#'
#' Curvilinear Component Analysis (CRCA) is a type of self-organizing algorithms for
#' manifold learning. Like MDS, it aims at minimizing a cost function (\emph{Stress})
#' based on pairwise proximity. Parameter \code{lambda} is a heaviside function for
#' penalizing distance pair of embedded data, and \code{alpha} controls learning rate
#' similar to that of subgradient method in that at each iteration \eqn{t} the gradient is
#' weighted by \eqn{\alpha /t}.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param lambda threshold value.
#' @param alpha initial value for updating.
#' @param maxiter maximum number of iterations allowed.
#' @param tolerance stopping criterion for maximum absolute discrepancy between two distance matrices.
#'
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{niter}{the number of iterations until convergence.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#'
#' @examples
#' \donttest{
#' ## generate sample data
#' X <- aux.gensamples(n=200)
#'
#' ## different initial learning rates
#' out1 <- do.crca(X,alpha=1)
#' out2 <- do.crca(X,alpha=5)
#' out3 <- do.crca(X,alpha=10)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="alpha=1.0")
#' plot(out2$Y, main="alpha=5.0")
#' plot(out3$Y, main="alpha=10.0")
#' par(opar)
#' }
#'
#'
#' @references
#' \insertRef{demartines_curvilinear_1997}{Rdimtools}
#'
#' \insertRef{herault_curvilinear_1999}{Rdimtools}
#'
#' @seealso \code{\link{do.crda}}
#' @author Kisung You
#' @rdname nonlinear_CRCA
#' @concept nonlinear_methods 
#' @export
do.crca <- function(X,ndim=2,lambda=1.0,alpha=1.0,maxiter=1000,tolerance=1e-6){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  vlnum = 1e+10
  # 1. data X
  aux.typecheck(X)
  # 2. ndim
  ndim = as.integer(ndim)
  if(!check_NumMM(ndim,1,ncol(X),compact=FALSE)){stop("* do.crca : 'ndim' is a positive integer in (1,#(covariates)).")}
  # 3. lambda
  lambda = as.double(lambda)
  if(!check_NumMM(lambda,0,vlnum,compact=FALSE)){stop("* do.crca : 'lambda' should be a positive real number.")}
  # 4. alpha
  alpha = as.double(alpha)
  if(!check_NumMM(alpha,0,vlnum,compact=FALSE)){stop("* do.crca : 'alpha' should be a positive real number.")}
  # 5. maxiter
  maxiter = as.integer(maxiter)
  if(!check_NumMM(maxiter,1,vlnum,compact=FALSE)){stop("* do.crca : 'maxiter' should be a not-so-small positive integer.")}
  # 6. tolerance
  tolerance = as.double(tolerance)
  if(!check_NumMM(alpha,0,vlnum,compact=FALSE)){stop("* do.crca : 'tolerance' should be a positive real number.")}

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. Pairwise Distance Matrix
  Xij = as.matrix(dist(X))
  #   2. Initialization via PCA
  Yinit = do.pca(X, ndim=ndim)$Y
  #   3. vecselector for random-number generation
  vecselector = as.vector(sample(0:(nrow(X)-1), maxiter, replace=TRUE))
  #   4. main computation
  Youtput = method_crca(Xij,Yinit,lambda,alpha,maxiter,tolerance,vecselector)

  #------------------------------------------------------------------------
  ## RETURN OUTPUT
  trfinfo = list()
  trfinfo$type = "null"
  trfinfo$algtype = "nonlinear"
  result = list()
  result$Y = Youtput$Y
  result$niter = Youtput$niter
  result$trfinfo = trfinfo
  return(result)
}
