#' Diversity-Induced Self-Representation
#'
#' Diversity-Induced Self-Representation (DISR) is a feature selection method that aims at
#' ranking features by both representativeness and diversity. Self-representation controlled by
#' \code{lbd1} lets the most representative features to be selected, while \code{lbd2} penalizes
#' the degree of inter-feature similarity to enhance diversity from the chosen features.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param lbd1 nonnegative number to control the degree of self-representation.
#' @param lbd2 nonnegative number to control the degree of feature similarity.
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
#' #### try different lbd combinations
#' out1 = do.disr(X, lbd1=1, lbd2=1)
#' out2 = do.disr(X, lbd1=1, lbd2=5)
#' out3 = do.disr(X, lbd1=5, lbd2=1)
#' out4 = do.disr(X, lbd1=5, lbd2=5)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2))
#' plot(out1$Y, main="(lbd1,lbd2)=(1,1)", col=label, pch=19)
#' plot(out2$Y, main="(lbd1,lbd2)=(1,5)", col=label, pch=19)
#' plot(out3$Y, main="(lbd1,lbd2)=(5,1)", col=label, pch=19)
#' plot(out4$Y, main="(lbd1,lbd2)=(5,5)", col=label, pch=19)
#' par(opar)
#' }
#'
#' @references
#' \insertRef{liu_unsupervised_2017}{Rdimtools}
#'
#' @seealso \code{\link{do.rsr}}
#' @author Kisung You
#' @rdname linear_DISR
#' @concept linear_methods
#' @export
do.disr <- function(X, ndim=2, preprocess=c("null","center","scale","cscale","whiten","decorrelate"),
                    lbd1=1.0, lbd2=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.disr : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. lbd1 and lbd2
  if ((length(lbd1)>1)||(lbd1<0)){
    stop("* do.disr : 'lbd1' should be a nonnegative real number.")
  }
  if ((length(lbd2)>1)||(lbd2<0)){
    stop("* do.disr : 'lbd2' should be a nonnegative real number.")
  }

  #------------------------------------------------------------------------
  ## COMPUTATION
  #   1. data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. main computation
  score      = method_disr(pX, lbd1, lbd2)
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
