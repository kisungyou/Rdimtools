#' Sparse Principal Component Analysis
#'
#' Sparse PCA (\code{do.spca}) is a variant of PCA in that each loading - or, principal
#' component - should be sparse. Unlike \code{\link[elasticnet]{spca}} from \pkg{elasticnet} package,
#' we only support sparsity type of \code{varnum} because it is more intuitive in explaining
#' up to how many elements can a loading have non-zero elements. We adopted semidefinite relaxation of
#' the original problem and optimization is done via \pkg{CVXR}.
#'
#' @param X an \code{(n-by-p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center" and "decorrelate", or "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param suppk the number of nonzero elements at each principal component.
#' Either a number or vector of length \code{ndim} can be supplied. Note that
#' every element in \code{suppk} should have range of [1,p].
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \code{(n-by-ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \code{(p-by-ndim)} whose columns are principal components.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' ## generate default dataset and make its dimension two-folds.
#' Xpart <- aux.gensamples()
#' X <- cbind(Xpart,Xpart)
#'
#' ## two different choice of 'suppk'
#' out1 <- do.spca(X,ndim=2,suppk=1)
#' out2 <- do.spca(X,ndim=2,suppk=6)
#'
#' ## Visualize principal components
#' par(mfrow=c(1,2))
#' image(out1$projection, main="support = 1")
#' image(out2$projection, main="support = 6")
#'
#' @references
#' \insertRef{zou_sparse_2006}{Rdimtools}
#'
#' \insertRef{daspremont_direct_2007}{Rdimtools}
#'
#' @seealso \code{\link{do.pca}}
#' @author Kisung You
#' @rdname linear_SPCA
#' @export
do.spca <- function(X,ndim=2,preprocess=c("center","decorrelate","whiten"),suppk=ceiling(ncol(X)/10)){
  ## PREPROCESSING
  # 1. data typecheck
  aux.typecheck(X)
  # 2. ndim
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.spca : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  # 3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  trfinfo$algtype = "linear"
  pX      = tmplist$pX
  # 4. suppk
  suppk = as.vector(suppk)
  if (length(round(suppk))==1){
    veck = rep(suppk,ndim)
  } else if (length(suppk)!=ndim){
    stop("* do.spca : 'suppk' as a vector input should have length of 'ndim'.")
  } else {
    veck = round(suppk)
  }
  if ((any(veck<1))||(any(veck>ncol(X)))){
    stop("* do.spca : 'suppk' should have values in [1,ncol(X)].")
  }

  ## MAIN COMPUTATION
  # 1. preliminary
  n = nrow(X) # obs : numbers
  p = ncol(X) # obs : dimension
  S = cov(pX) # sample covariance
  projection = array(0,c(p,ndim))

  # 2. iteration for CVXR
  for (i in 1:ndim){
    # 2-1. for each dimension, we have different support condition
    tgtk   = veck[i]
    # 2-2. CVXR specification
    V      = CVXR::Semidef(p)
    obj    = matrix_trace(S%*%V)
    constr = list(matrix_trace(V)==1,
                  sum(abs(V))<=tgtk
    )
    prob   = Problem(Maximize(obj), constr)
    solcvx = solve(prob)
    # 2-3. extract
    eigV = eigen(matrix(solcvx$getValue(V), nrow=p))
    v    = as.matrix(eigV$vectors[,1]*as.double(sqrt(eigV$values[1])))
    # 2-4. record in Y
    projection[,i] = as.vector(v)
    # 2-5. update S
    S = (as.double(t(v)%*%S%*%v))*outer(as.vector(v),as.vector(v))
  }

  ## Return output
  result   = list()
  result$Y = pX%*%projection
  result$projection = projection
  result$trfinfo = trfinfo
  return(result)
}
