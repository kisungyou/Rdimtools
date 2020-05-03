#' Non-convex Regularized Self-Representation
#'
#' In the standard, convex RSR problem (\code{\link{do.rsr}}), row-sparsity for self-representation is
#' acquired using matrix \eqn{\ell_{2,1}} norm, i.e, \eqn{\|W\|_{2,1} = \sum \|W_{i:}\|_2}. Its non-convex
#' extension aims at achieving higher-level of sparsity using arbitrarily chosen \eqn{\|W\|_{2,l}} norm for
#' \eqn{l\in (0,1)} and this exploits Iteratively Reweighted Least Squares (IRLS) algorithm for computation.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param expl an exponent in \eqn{\ell_{2,l}} norm for sparsity. Must be in \eqn{(0,1)}, or \eqn{l=1} reduces to RSR problem.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param lbd nonnegative number to control the degree of self-representation by imposing row-sparsity.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
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
#' #### try different exponents for regularization
#' out1 = do.nrsr(X, expl=0.01)
#' out2 = do.nrsr(X, expl=0.1)
#' out3 = do.nrsr(X, expl=0.5)
#'
#' #### visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=label, main="NRSR::expl=0.01")
#' plot(out2$Y, pch=19, col=label, main="NRSR::expl=0.1")
#' plot(out3$Y, pch=19, col=label, main="NRSR::expl=0.5")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zhu_nonconvex_2017}{Rdimtools}
#'
#' @seealso \code{\link{do.rsr}}
#' @author Kisung You
#' @rdname linear_NRSR
#' @concept linear_methods
#' @export
do.nrsr <- function(X, ndim=2, expl=0.5, preprocess=c("null","center","scale","cscale","whiten","decorrelate"), lbd=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.nrsr : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. lbd
  if ((length(lbd)>1)||(lbd<0)){
    stop("* do.nrsr : 'lbd' should be a nonnegative real number.")
  }
  #   5. expl
  if ((length(expl)>1)||(expl<=0)||(expl>1)){
    stop("* do.nrsr : 'expl' must be a real number in (0,1). 'expl=1' reduces to the convex 'do.rsr' problem.")
  }

  #------------------------------------------------------------------------
  ## COMPUTATION
  #   1. data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. main computation
  score      = method_nrsr(pX, lbd, sqrt(10*.Machine$double.eps), expl)
  idxvec     = base::order(score, decreasing=TRUE)[1:ndim]
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$featidx = idxvec
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
