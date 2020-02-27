#' Enhanced Locality Preserving Projection (2013)
#'
#' Enhanced Locality Preserving Projection proposed in 2013 (ELPP2) is built upon
#' a parameter-free philosophy from PFLPP. It further aims to exclude its projection
#' to be uncorrelated in the sense that the scatter matrix is placed in a generalized eigenvalue problem.
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
#' \dontrun{
#' ## generate swiss roll data
#' X = aux.gensamples(n=200)
#'
#' ## compare with PCA and PFLPP
#' out1 = do.pca(X, ndim=2)
#' out2 = do.pflpp(X, ndim=2)
#' out3 = do.elpp2(X, ndim=2)
#'
#' ## visualize
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(out1$Y, main="PCA")
#' plot(out2$Y, main="Parameter-Free LPP")
#' plot(out3$Y, main="Enhanced LPP (2013)")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{dornaika_enhanced_2013}{Rdimtools}
#'
#' @seealso \code{\link{do.pflpp}}
#' @author Kisung You
#' @rdname linear_ELPP2
#' @export
do.elpp2 <- function(X, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten")){
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
  #      unlike PFLPP, now we use scatter matrix for RHS
  LHS     = t(pX)%*%L%*%pX
  meanvec = as.vector(colMeans(pX))
  S       = aux_scatter(pX, meanvec)

  #   4. projection vectors
  projection = aux.geigen(LHS, S, ndim, maximal=FALSE)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
