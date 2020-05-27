#' Locality and Similarity Preserving Embedding
#'
#' Locality and Similarity Preserving Embedding (LSPE) is a feature selection method based on Neighborhood Preserving Embedding (\code{\link{do.npe}}) and
#' Sparsity Preserving Projection (\code{\link{do.spp}}) by first building a neighborhood graph and
#' then mapping the locality structure to reconstruct coefficients such that data similarity is preserved.
#' Use of \eqn{\ell_{2,1}} norm boosts to impose column-sparsity that enables feature selection procedure.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param alpha nonnegative number to control \eqn{\ell_{2,1}} norm of projection.
#' @param beta nonnegative number to control the degree of local similarity.
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
#' \donttest{
#' #### generate R12in72 dataset
#' set.seed(100)
#' X = aux.gensamples(n=50, dname="R12in72")
#'
#' #### try different bandwidth values
#' out1 = do.lspe(X, bandwidth=0.1)
#' out2 = do.lspe(X, bandwidth=1)
#' out3 = do.lspe(X, bandwidth=10)
#'
#' #### visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="LSPE::bandwidth=0.1")
#' plot(out2$Y, main="LSPE::bandwidth=1")
#' plot(out3$Y, main="LSPE::bandwidth=10")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{fang_locality_2014}{Rdimtools}
#'
#' @seealso \code{\link{do.rsr}}
#' @rdname feature_LSPE
#' @author Kisung You
#' @concept feature_methods
#' @export
do.lspe <- function(X, ndim=2, preprocess=c("null","center","scale","cscale","whiten","decorrelate"),
                    alpha=1.0, beta=1.0, bandwidth=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lspe : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. alpha and beta
  if ((length(alpha)>1)||(alpha<0)){
    stop("* do.lspe : 'alpha' should be a nonnegative real number.")
  }
  if ((length(beta)>1)||(beta<0)){
    stop("* do.lspe : 'beta' should be a nonnegative real number.")
  }
  if ((length(bandwidth)>1)||(bandwidth<=0)){
    stop("* do.lspe : 'bandwidth' should be a nonnegative real number.")
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
  L       = D-S
  #   3. main part
  tpX     = t(pX)
  score   = method_lspe(tpX, ndim, alpha, beta, L)
  idxvec  = base::order(score, decreasing=TRUE)[1:ndim]
  #   4. get the projection based on top scored ones
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
