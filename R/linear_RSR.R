#' Regularized Self-Representation
#'
#' Given a data matrix \eqn{X} where observations are stacked in a row-wise manner,
#' Regularized Self-Representation (RSR) aims at finding a solution to following optimization problem
#' \deqn{\textrm{min}~ \|X-XW\|_{2,1} + \lambda \| W \|_{2,1}}
#' where \eqn{\|W\|_{2,1} = \sum_{i=1}^{m} \|W_{i:} \|_2} is an \eqn{\ell_{2,1}} norm that imposes
#' row-wise sparsity constraint.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
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
#' X     = as.matrix(iris[,1:4])
#' label = as.integer(iris$Species)
#'
#' #### try different lbd combinations
#' out1 = do.rsr(X, lbd=0.1)
#' out2 = do.rsr(X, lbd=1)
#' out3 = do.rsr(X, lbd=10)
#'
#' #### visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(3,1))
#' plot(out1$Y, col=label, main="RSR::lbd=0.1")
#' plot(out2$Y, col=label, main="RSR::lbd=1")
#' plot(out3$Y, col=label, main="RSR::lbd=10")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zhu_unsupervised_2015}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_RSR
#' @export
do.rsr <- function(X, ndim=2, preprocess=c("null","center","scale","cscale","whiten","decorrelate"), lbd=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.rsr : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. lbd
  if ((length(lbd)>1)||(lbd<0)){
    stop("* do.rsr : 'lbd' should be a nonnegative real number.")
  }

  #------------------------------------------------------------------------
  ## COMPUTATION
  #   1. data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. main computation
  score      = method_rsr(pX, lbd, sqrt(10*.Machine$double.eps))
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
