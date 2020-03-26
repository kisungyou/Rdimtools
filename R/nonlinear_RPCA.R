#' Robust Principal Component Analysis
#'
#' Robust PCA (RPCA) is not like other methods in this package as finding explicit low-dimensional embedding with reduced number of columns.
#' Rather, it is more of a decomposition method of data matrix \eqn{X}, possibly noisy, into low-rank and sparse matrices by
#' solving the following,
#' \deqn{\textrm{minimize}\quad \|L\|_* + \lambda \|S\|_1 \quad{s.t.} L+S=X}
#' where \eqn{L} is a low-rank matrix, \eqn{S} is a sparse matrix and \eqn{\|\cdot\|_*} denotes nuclear norm, i.e., sum of singular values. Therefore,
#' it should be considered as \emph{preprocessing} procedure of denoising. Note that after RPCA is applied, \eqn{L} should be used
#' as kind of a new data matrix for any manifold learning scheme to be applied.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param mu an augmented Lagrangian parameter
#' @param lambda parameter for the sparsity term \eqn{\|S\|_1}. Default value is given accordingly to the referred paper.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{L}{an \eqn{(n\times p)} low-rank matrix.}
#' \item{S}{an \eqn{(n\times p)} sparse matrix.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' ## Load Iris data and put some noise
#' data(iris)
#' noise = 0.2
#' X = as.matrix(iris[,1:4])
#' X = X + matrix(noise*rnorm(length(X)), nrow=nrow(X))
#' label = as.integer(iris$Species)
#'
#' ## try different regularization parameters
#' rpca1 = do.rpca(X, lambda=0.1)
#' rpca2 = do.rpca(X, lambda=1)
#' rpca3 = do.rpca(X, lambda=10)
#'
#' ## apply identical PCA methods
#' out1 = do.pca(rpca1$L, ndim=2)
#' out2 = do.pca(rpca2$L, ndim=2)
#' out3 = do.pca(rpca3$L, ndim=2)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, main="RPCA+PCA::lambda=0.1")
#' plot(out2$Y, col=label, main="RPCA+PCA::lambda=1")
#' plot(out3$Y, col=label, main="RPCA+PCA::lambda=10")
#' par(opar)
#'
#' @references
#' \insertRef{candes_robust_2011}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_RPCA
#' @export
do.rpca <- function(X, mu=1.0, lambda=sqrt(1/(max(dim(X)))),
                    preprocess=c("null","center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. mu and lambda
  muval  = as.double(mu)
  lbdval = as.double(lambda)
  if (!check_NumMM(muval,0,Inf,compact=FALSE)){stop("* do.rpca : 'mu' should be a POSITIVE real number.")}
  if (!check_NumMM(lbdval,0,Inf,compact=FALSE)){stop("* do.rpca : 'lambda' should be a POSITIVE real number.")}

  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION USING ADMM PACKAGE
  admmrun = ADMM::admm.rpca(pX,lambda=lbdval,mu=muval)

  #------------------------------------------------------------------------
  ## REPORT RESULTS
  result = list()
  result$L = admmrun$L
  result$S = admmrun$S
  result$trfinfo = trfinfo
  return(result)
}
