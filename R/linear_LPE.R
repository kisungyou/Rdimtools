#' Locality Pursuit Embedding
#'
#' Locality Pursuit Embedding (LPE) is an unsupervised linear dimension reduction method.
#' It aims at preserving local structure by solving a variational problem that models
#' the local geometrical structure by the Euclidean distances.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center" and other options of "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param numk size of \eqn{k}-nn neighborhood in original dimensional space.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate swiss roll with auxiliary dimensions
#' n = 100
#' theta = runif(n)
#' h     = runif(n)
#' t     = (1+2*theta)*(3*pi/2)
#' X     = array(0,c(n,10))
#' X[,1] = t*cos(t)
#' X[,2] = 21*h
#' X[,3] = t*sin(t)
#' X[,4:10] = matrix(runif(7*n), nrow=n)
#'
#' ## try with different neighborhood sizes
#' out1 = do.lpe(X, numk=5)
#' out2 = do.lpe(X, numk=10)
#' out3 = do.lpe(X, numk=25)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="LPE::numk=5")
#' plot(out2$Y[,1], out2$Y[,2], main="LPE::numk=10")
#' plot(out3$Y[,1], out3$Y[,2], main="LPE::numk=25")
#' }
#'
#' @references
#' \insertRef{min_locality_2004}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LPE
#' @export
do.lpe <- function(X, ndim=2, preprocess=c("center","decorrelate","whiten"), numk=max(ceiling(nrow(X)/10),2)){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lpe : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. numk
  numk = as.integer(numk)
  if (!check_NumMM(numk,1,n/2,compact=FALSE)){stop("* do.lpe : 'numk' should be an integer in [2,nrow(X)/2).")}
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
  #   2. neighborhood creation
  nbdtype = c("knn",numk)
  nbdsymmetric = "asymmetric"
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR LPE
  #   1. build L
  L = array(0,c(n,n))
  onesN = array(1,c(n,n))
  for (i in 1:n){
    vecdi = (as.vector(nbdmask[i,])*1.0)
    K     = sum(vecdi)
    Di    = diag(vecdi)
    L     = L + Di + ((1/K)*(Di%*%onesN%*%Di))
  }
  #   2. find cost function
  costTop = t(pX)%*%L%*%pX
  #   3. find projection matrix
  projection = aux.adjprojection(RSpectra::eigs(costTop, ndim)$vectors)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
