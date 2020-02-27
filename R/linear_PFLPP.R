#' Parameter-Free Locality Preserving Projection
#'
#' Conventional LPP is known to suffer from sensitivity upon choice of parameters, especially
#' in building neighborhood information. Parameter-Free LPP (PFLPP) takes an alternative step
#' to use normalized Pearson correlation, taking an average of such similarity as a threshold
#' to decide which points are neighbors of a given datum.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' ## generate swiss roll data
#' X = aux.gensamples(n=200)
#'
#' ## compare with PCA
#' out1 = do.pca(X, ndim=2)
#' out2 = do.pflpp(X, ndim=2)
#'
#' ## visualize
#' opar <- par(mfrow=c(1,2), no.readonly=TRUE)
#' plot(out1$Y, main="PCA")
#' plot(out2$Y, main="Parameter-Free LPP")
#' par(opar)
#'
#' @references
#' \insertRef{dornaika_enhanced_2013}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_PFLPP
#' @export
do.pflpp <- function(X, ndim=2, preprocess=c("center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.pflpp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. pearson correlation among samples
  tmpP = stats::cor(t(pX)); mintmpP = as.double(min(tmpP))

  #   3. normalized correlation
  P = (tmpP-mintmpP)/(1-mintmpP)

  #   4. mean for each datum
  vecm = rowMeans(P)

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR PFLPP
  #   1. build adjacency
  A = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (P[i,j] > max(vecm[i],vecm[j])){
        A[i,j] = P[i,j]
        A[j,i] = P[i,j]
      }
    }
  }
  #   2. build for L
  D = diag(rowSums(A))
  L = D-A
  #   3. cost function with respect to geigen problem
  LHS = t(pX)%*%L%*%pX
  RHS = t(pX)%*%D%*%pX
  #   4. projection vectors
  projection = aux.geigen(LHS, RHS, ndim, maximal=FALSE)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
