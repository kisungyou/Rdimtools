#' Adaptive Subspace Iteration
#'
#' Adaptive Subspace Iteration (ASI) iteratively finds the best subspace to perform data clustering. It can be regarded as
#' one of remedies for clustering in high dimensional space. Eigenvectors of a within-cluster scatter matrix are used
#' as basis of projection.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
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
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150, 50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' ## try different tolerance level
#' out1 = do.asi(X, abstol=1e-2)
#' out2 = do.asi(X, abstol=1e-3)
#' out3 = do.asi(X, abstol=1e-4)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="ASI::tol=1e-2", pch=19, col=label)
#' plot(out2$Y, main="ASI::tol=1e-3", pch=19, col=label)
#' plot(out3$Y, main="ASI::tol=1e-4", pch=19, col=label)
#' par(opar)
#'
#' @references
#' \insertRef{li_document_2004}{Rdimtools}
#'
#' @seealso  \code{\link{do.ldakm}}
#' @author Kisung You
#' @rdname linear_ASI
#' @concept linear_methods
#' @export
do.asi <- function(X, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"),
                   maxiter=10, abstol=1e-3){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim as 'd' and 'k' the number of clusters
  d = as.integer(ndim)
  if (!check_ndim(d,p)){stop("* do.asi : 'ndim' is a positive integer in [1,#(covariates)).")}
  k = as.integer(d+1)
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. maxiter
  maxiter = as.integer(maxiter)
  if (!check_NumMM(maxiter,3,1000000)){stop("* do.asi : 'maxiter' should be a large positive integer.")}
  #   5. abstol
  abstol = as.double(abstol)
  if (!check_NumMM(abstol,0,0.5,compact=FALSE)){stop("* do.asi : 'abstol' should be a small nonnegative number for stopping criterion.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
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
    # Sb = M%*%t(H)%*%H%*%t(M)
    # 3-3. BRANCHING :: Solve for Eigenvectors
    Unew = aux.adjprojection(RSpectra::eigs(Sw,ndim,which="SR")$vectors)
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
