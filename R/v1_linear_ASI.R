#' Adaptive Subspace Iteration
#'
#' Adaptive Subspace Iteration (ASI) iteratively finds the best subspace to perform data clustering. It can be regarded as
#' one of remedies for clustering in high dimensional space. Eigenvectors of a within-cluster scatter matrix are used
#' as basis of projection.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations (default: 100).}
#' \item{abstol}{absolute tolerance stopping criterion (default: 1e-8).}
#' }
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
#' ## use iris data
#' data(iris, package="Rdimtools")
#' set.seed(100)
#' subid = sample(1:150, 50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' ## compare ASI with other methods
#' outASI = do.asi(X)
#' outPCA = do.pca(X)
#' outLDA = do.lda(X, label)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(outASI$Y, pch=19, col=label, main="ASI")
#' plot(outPCA$Y, pch=19, col=label, main="PCA")
#' plot(outLDA$Y, pch=19, col=label, main="LDA")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{li_document_2004}{Rdimtools}
#'
#' @seealso  \code{\link{do.ldakm}}
#' @author Kisung You
#' @rdname linear_ASI
#' @concept linear_methods
#' @export
do.asi <- function(X, ndim=2, ...){
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
  # #   3. preprocess
  # if (missing(preprocess)){
  #   algpreprocess = "center"
  # } else {
  #   algpreprocess = match.arg(preprocess)
  # }
  # #   4. maxiter
  # maxiter = as.integer(maxiter)
  # if (!check_NumMM(maxiter,3,1000000)){stop("* do.asi : 'maxiter' should be a large positive integer.")}
  # #   5. abstol
  # abstol = as.double(abstol)
  # if (!check_NumMM(abstol,0,0.5,compact=FALSE)){stop("* do.asi : 'abstol' should be a small nonnegative number for stopping criterion.")}

  # Extra parameters
  params  = list(...)
  pnames  = names(params)

  if ("abstol"%in%pnames){
    abstol = max(.Machine$double.eps, as.double(params$abstol))
  } else {
    abstol = 10^(-8)
  }
  if ("maxiter"%in%pnames){
    maxiter = max(5, round(params$maxiter))
  } else {
    maxiter = 100
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data
  # tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  # trfinfo = tmplist$info
  # pX      = tmplist$pX
  #   2. initialize
  Uold = ldakm_PCAbasis(X, ndim)
  #   3. iterate
  incstop = 10.0
  citer   = 1
  while (incstop > abstol){
    # 3-1. LDA-KM(1) : k-means in projected space
    projected = X%*%Uold
    pXkmeans  = kmeans(projected, k)
    # 3-2. LDA-KM(2) : learn again
    # 1. build H
    H = ldakm_BuildH(pXkmeans$cluster)     # H : (n-times-k)
    M = (t(X)%*%H%*%aux.pinv(t(H)%*%H))   # M : (p-times-k)
    # 2. build Sw (p-by-p)
    Swterm1 = t(X)-(M%*%t(H))
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
  result$Y = X%*%projection
  result$projection = projection
  result$algorithm  = "linear:ASI"
  return(structure(result, class="Rdimtools"))
}
