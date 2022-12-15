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
#' @param lbd nonnegative number to control the degree of self-representation by imposing row-sparsity.
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' \donttest{
#' ## load iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' #### try different lbd combinations
#' out1 = do.rsr(X, lbd=0.1)
#' out2 = do.rsr(X, lbd=1)
#' out3 = do.rsr(X, lbd=10)
#'
#' #### visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=label, main="RSR::lbd=0.1")
#' plot(out2$Y, pch=19, col=label, main="RSR::lbd=1")
#' plot(out3$Y, pch=19, col=label, main="RSR::lbd=10")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zhu_unsupervised_2015}{Rdimtools}
#'
#' @author Kisung You
#' @rdname feature_RSR
#' @concept feature_methods
#' @export
do.rsr <- function(X, ndim=2, lbd=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.rsr : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. lbd
  if ((length(lbd)>1)||(lbd<0)){
    stop("* do.rsr : 'lbd' should be a nonnegative real number.")
  }

  #------------------------------------------------------------------------
  ## COMPUTATION
  score      = method_rsr(X, lbd, sqrt(10*.Machine$double.eps))
  idxvec     = base::order(score, decreasing=TRUE)[1:ndim]
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = X%*%projection
  result$featidx = idxvec
  result$projection = projection
  result$algorithm = "linear:RSR"
  return(structure(result, class="Rdimtools"))
}
