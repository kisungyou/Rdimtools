#' Structure Preserving Unsupervised Feature Selection
#'
#' This unsupervised feature selection method is based on self-expression model, which means that
#' the cost function involves difference in self-representation. It does not explicitly
#' require learning the clusterings and different features are weighted individually
#' based on their relative importance. The cost function involves two penalties,
#' sparsity and preservation of local structure.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param alpha nonnegative number to control sparsity in rows of matrix of representation coefficients.
#' @param beta nonnegative number to control the degree of local-structure preservation.
#' @param bandwidth positive number for Gaussian kernel bandwidth to define similarity.
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
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.integer(iris$Species)
#'
#' #### try different bandwidth values
#' out1 = do.spufs(X, bandwidth=0.1)
#' out2 = do.spufs(X, bandwidth=1)
#' out3 = do.spufs(X, bandwidth=10)
#'
#' #### visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, main="SPUFS::bandwidth=0.1")
#' plot(out2$Y, col=label, main="SPUFS::bandwidth=1")
#' plot(out3$Y, col=label, main="SPUFS::bandwidth=10")
#' par(opar)
#'
#' @references
#' \insertRef{lu_structure_2018}{Rdimtools}
#'
#' @rdname linear_SPUFS
#' @author Kisung You
#' @concept linear_methods 
#' @export
do.spufs <- function(X, ndim=2, preprocess=c("null","center","scale","cscale","whiten","decorrelate"),
                     alpha=1.0, beta=1.0, bandwidth=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.spufs : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. alpha and beta
  if ((length(alpha)>1)||(alpha<0)){
    stop("* do.spufs : 'alpha' should be a nonnegative real number.")
  }
  if ((length(beta)>1)||(beta<0)){
    stop("* do.spufs : 'beta' should be a nonnegative real number.")
  }
  if ((length(bandwidth)>1)||(bandwidth<=0)){
    stop("* do.spufs : 'bandwidth' should be a nonnegative real number.")
  }

  #------------------------------------------------------------------------
  ## COMPUTATION
  #   1. data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #   2. weight matrix and graph laplacian
  Dmat    = as.matrix(dist(pX)); colnames(Dmat)=NULL; rownames(Dmat)=NULL;
  S       = exp(-(Dmat^2)/bandwidth)
  diag(S) = 0 # (n x n) matrix
  D       = diag(rowSums(S))
  Ls      = D-S
  #   3. weight-vector score and find idx
  epsilon = sqrt(10*.Machine$double.eps)
  score   = method_spufs(pX, Ls, alpha, beta, epsilon)
  idxvec  = base::order(score, decreasing=TRUE)[1:ndim]
  #   4. find projection
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
