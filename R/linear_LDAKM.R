#' Combination of LDA and K-means
#'
#' \code{do.ldakm} is an unsupervised subspace discovery method that combines linear discriminant analysis (LDA) and K-means algorithm.
#' It tries to build an adaptive framework that selects the most discriminative subspace. It iteratively applies two methods in that
#' the clustering process is integrated with the subspace selection, and continuously updates its discrimative basis. From its formulation
#' with respect to generalized eigenvalue problem, it can be considered as generalization of Adaptive Subspace Iteration (ASI) and Adaptive Dimension Reduction (ADR).
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center" and other options of "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param maxiter maximum number of iterations allowed.
#' @param abstol stopping criterion for incremental change in projection matrix.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## generate swiss roll data
#' X = aux.gensamples(n=123)
#'
#' ## try different tolerance level
#' out1 = do.ldakm(X, abstol=1e-2)
#' out2 = do.ldakm(X, abstol=1e-3)
#' out3 = do.ldakm(X, abstol=1e-4)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="LDA-KM::tol=1e-2")
#' plot(out2$Y[,1], out2$Y[,2], main="LDA-KM::tol=1e-3")
#' plot(out3$Y[,1], out3$Y[,2], main="LDA-KM::tol=1e-4")
#'
#' @references
#' \insertRef{ding_adaptive_2007}{Rdimtools}
#' @seealso \code{\link{do.asi}}, \code{\link{do.adr}}
#' @author Kisung You
#' @rdname linear_LDAKM
#' @export
do.ldakm <- function(X, ndim=2, preprocess=c("center","decorrelate","whiten"), maxiter=10, abstol=1e-3){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim as 'd' and 'k' the number of clusters
  d = as.integer(ndim)
  if (!check_ndim(d,p)){stop("* do.ldakm : 'ndim' is a positive integer in [1,#(covariates)).")}
  k = as.integer(d+1)
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. maxiter
  maxiter = as.integer(maxiter)
  if (!check_NumMM(maxiter,3,1000000)){stop("* do.ldakm : 'maxiter' should be a large positive integer.")}
  #   5. abstol
  abstol = as.double(abstol)
  if (!check_NumMM(abstol,0,0.5,compact=FALSE)){stop("* do.ldakm : 'abstol' should be a small nonnegative number for stopping criterion.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
  #   2. initialize
  Uold = ldakm_PCAbasis(pX, ndim)
  #   3. iterate
  incstop = 10.0
  citer   = 1
  while (incstop > abstol){
    # 3-1. LDA-KM(1) : k-means in projected space
    projected = pX%*%Uold
    pXkmeans  = kmeans(projected, k)
    # 3-2. LDA-KM(2) : learn again
    # 1. build H
    H = ldakm_BuildH(pXkmeans$cluster)     # H : (n-times-k)
    M = (t(pX)%*%H%*%aux.pinv(t(H)%*%H))   # M : (p-times-k)
    # 2. build Sw (p-by-p)
    Swterm1 = t(pX)-(M%*%t(H))
    Sw = Swterm1%*%t(Swterm1)
    # 3. build Sb (p-by-p)
    Sb = M%*%t(H)%*%H%*%t(M)
    # 3-3. BRANCHING :: Solve for Eigenvectors
    Unew = aux.geigen(Sb, Sw, ndim, maximal=TRUE)
    # 3-4. update
    incstop = base::norm(Uold-Unew,"f")
    citer = citer + 1
    Uold  = Unew
    if (citer >= maxiter){
      break
    }
  }
  #   4. we finally have projection
  projection = aux.adjprojection(Uold)


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}



#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
ldakm_PCAbasis <- function(X, ndim){
  basis = aux.adjprojection(RSpectra::eigs(cov(X), ndim)$vectors)
  return(basis)
}
# Build H = n-times-k matrix
#' @keywords internal
#' @noRd
ldakm_BuildH <- function(labeling){
  k = length(unique(labeling))
  n = length(labeling)
  H = array(0,c(n,k))
  for (i in 1:n){
    H[i,as.integer(labeling[i])] = 1
  }
  return(H)
}
