#' Linear Discriminant Analysis
#'
#' Linear Discriminant Analysis (LDA) originally aims to find a set of features
#' that best separate groups of data. Since we need \emph{label} information,
#' LDA belongs to a class of supervised methods of performing classification.
#' However, since it is based on finding \emph{suitable} projections, it can still
#' be used to do dimension reduction. We support both binary and multiple-class cases.
#' Note that the target dimension \code{ndim} should be \emph{less than or equal to} \code{K-1},
#' where \code{K} is the number of classes, or \code{K=length(unique(label))}. Our code
#' automatically gives bounds on user's choice to correspond to what theory has shown. See
#' the comments section for more details.
#'
#' @section Limit of Target Dimension Selection:
#' In unsupervised algorithms, selection of \code{ndim} is arbitrary as long as
#' the target dimension is lower-dimensional than original data dimension, i.e., \code{ndim < p}.
#' In LDA, it is \emph{not allowed}. Suppose we have \code{K} classes, then its formulation on
#' \eqn{S_B}, between-group variance, has maximum rank of \code{K-1}. Therefore, the maximal
#' subspace can only be spanned by at most \code{K-1} orthogonal vectors.
#'
#' @param X an \eqn{(n\times p)} matrix whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris dataset
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' lab   = as.factor(iris[,5])
#'
#' ## compare with PCA
#' outLDA = do.lda(X, lab, ndim=2)
#' outPCA = do.pca(X, ndim=2)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(outLDA$Y, col=lab, pch=19, main="LDA")
#' plot(outPCA$Y, col=lab, pch=19, main="PCA")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{fisher_use_1936}{Rdimtools}
#'
#' \insertRef{fukunaga_introduction_1990}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LDA
#' @concept linear_methods
#' @export
do.lda <- function(X, label, ndim=2){
  #------------------------------------------------------------------------
  ## BASIC
  #  data
  if (!is.matrix(X)){
    stop("* do.lda : 'X' should be a matrix.")
  }
  n = nrow(X)
  p = ncol(X)

  #  label
  label  = check_label(label, n)
  ulabel = unique(label)
  K      = length(ulabel)
  if (K==1){
    stop("* do.lda : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    warning("* do.lda : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }

  #  ndim
  if (!check_ndim(ndim,p)){
    stop("* do.lda : 'ndim' should be a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  if (ndim>=K){
    warning("* do.lda : by the nature of LDA, target dimension 'ndim' needs to be adjusted to match maximally permissible subspace.")
  }

  #------------------------------------------------------------------------
  ## Rewriting LDA
  # 1. split data into the list
  datlist = list()
  for (i in 1:length(ulabel)){
    datlist[[i]] = X[which(label==ulabel[i]),]
  }
  # 2. compute two types of scatter matrices
  #   2-1. error matrix/E/error variance
  scattermat <- function(x){
    return(cov(x)*(nrow(x)-1))
  }
  matE = array(0,c(p,p))
  for (i in 1:length(ulabel)){
    matE = matE + cov(datlist[[i]])*(nrow(datlist[[i]])-1)
  }
  matH = array(0,c(p,p))
  meanlist = lapply(datlist, colMeans)
  meantott = colMeans(X)
  for (i in 1:length(ulabel)){
    meandiff = as.vector(meanlist[[i]]-meantott)
    matH = matH + nrow(datlist[[i]])*outer(meandiff,meandiff)
  }

  W = aux.traceratio(matH, matE, ndim, 1e-6, 123)

  #------------------------------------------------------------------------
  # return results
  result = list()
  result$Y = X%*%W
  result$projection = W
  result$algorithm  = "linear:lda"
  return(structure(result, class="Rdimtools"))
}

#'
#' ## sum of outer products to calculate S_B and S_W
#' #' @keywords internal
#' #' @noRd
#' lda_outer <- function(X){
#'   p = ncol(X)
#'   output = array(0,c(p,p))
#'   for (i in 1:nrow(X)){
#'     output = output + outer(X[i,],X[i,])
#'   }
#'   return(output)
#' }
