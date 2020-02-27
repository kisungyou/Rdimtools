#' Robust Principal Component Analysis via Geometric Median
#'
#' This function robustifies the traditional PCA via an idea of geometric median.
#' To describe, the given data is first split into \code{k} subsets for each sample
#' covariance is attained. According to the paper, the median covariance is computed
#' under Frobenius norm and projection is extracted from the largest eigenvectors.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param k the number of subsets for \code{X} to be divided.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
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
#' ## try different numbers for subsets
#' out1 = do.rpcag(X, ndim=2, k=2)
#' out2 = do.rpcag(X, ndim=2, k=5)
#' out3 = do.rpcag(X, ndim=2, k=10)
#'
#' ## visualize
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(out1$Y, col=label, main="RPCAG::k=2")
#' plot(out2$Y, col=label, main="RPCAG::k=5")
#' plot(out3$Y, col=label, main="RPCAG::k=10")
#' par(opar)
#'
#' @references
#' \insertRef{minsker_geometric_2015}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_RPCAG
#' @export
do.rpcag <- function(X, ndim=2, k=5,
                     preprocess=c("center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.rpcag : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # 2. grouping
  k = round(k)
  n = nrow(X)
  if (n/k <= 2){
    warning("* do.rpcag : 'k' should be smaller than 'nrow(X)/2'.")
    k = max(min(k-1, floor(n/2)),2)
  }

  # 2. process : data preprocessing
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION
  projection = rpcag_internal(pX, k, ndim)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}


# sub-function ------------------------------------------------------------
#' @keywords internal
rpcag_internal <- function(X, k, ndim){
  n = base::nrow(X)
  p = base::ncol(X)

  # stack covariances
  covecs   = c()
  subsetid = aux.subsetid(n,k)
  for (i in 1:k){
    covecs = base::rbind(covecs,as.vector(stats::cov(X[subsetid[[i]],])))
  }

  # compute L1-median
  cov = matrix(as.vector(maotai::weiszfeld(covecs)), nrow=p)
  cov = (cov + t(cov))/2

  # compute top eigenvectors
  return(base::eigen(cov, symmetric = TRUE)$vectors[,1:ndim])
}


