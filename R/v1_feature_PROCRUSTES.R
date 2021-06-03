#' Feature Selection using PCA and Procrustes Analysis
#'
#' \code{do.procrustes} selects a set of features that best aligns PCA's coordinates in the embedded low dimension.
#' It iteratively selects each variable that minimizes Procrustes distance between configurations.
#'
#' @param X an \eqn{(n\times p)} matrix whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param intdim intrinsic dimension of PCA to be applied. It should be smaller than \code{ndim}.
#' @param cor mode of eigendecomposition. \code{FALSE} for decomposing covariance, and \code{TRUE} for correlation matrix in PCA.
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' ## it is known that feature 3 and 4 are more important.
#' data(iris)
#' iris.dat = as.matrix(iris[,1:4])
#' iris.lab = as.factor(iris[,5])
#'
#' ## try different strategy
#' out1 = do.procrustes(iris.dat, cor=TRUE)
#' out2 = do.procrustes(iris.dat, cor=FALSE)
#' out3 = do.mifs(iris.dat, iris.lab, beta=0)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1, 3))
#' plot(out1$Y, pch=19, col=iris.lab, main="PCA with Covariance")
#' plot(out2$Y, pch=19, col=iris.lab, main="PCA with Correlation")
#' plot(out3$Y, pch=19, col=iris.lab, main="MIFS")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{krzanowski_selection_1987a}{Rdimtools}
#'
#' @rdname feature_PROCRUSTES
#' @author Kisung You
#' @concept feature_methods
#' @export
do.procrustes <- function(X, ndim=2, intdim=(ndim-1), cor=TRUE){
  #------------------------------------------------------------------------
  ## Basic
  if (!is.matrix(X)){
    stop("* do.procrustes : 'X' should be a matrix.")
  }
  ndim = min(max(1, round(ndim)), ncol(X)-1)
  n    = nrow(X)
  p    = ncol(X)
  if (!check_ndim(ndim,p)){
    stop("* do.procrustes : 'ndim' is a positive integer in [1,#(covariates)].")
  }

  ## Others
  intdim = round(intdim)
  if (intdim >= ndim){
    stop("* do.procrustes : intrinsic dimension parameter 'intdim' should be smaller than 'dim'.")
  }
  mycor = as.logical(cor)

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  pY = as.matrix(dt_pca(X, intdim, mycor)$Y)

  #------------------------------------------------------------------------
  ## COMPUTATION : ITERATIVE SELECTION
  idall = 1:p # all indices
  idsel = c() # selected  ones
  idcan = 1:p # candidate ones

  for (it in 1:ndim){
    # candidate
    can.n    = length(idcan)
    can.cost = rep(0,can.n)
    for (i in 1:can.n){
      tilX = X[,setdiff(idall, c(idsel, idcan[i]))]
      tilZ = dt_pca(tilX, intdim, mycor)$Y
      tilZ = matrix(base::scale(tilZ, center=TRUE, scale=FALSE), ncol=intdim)
      can.cost[i] = procrustes_cost(pY, tilZ, intdim)
    }
    # the one with the minimal cost is the one to be chosen
    fopt = idcan[which.min(can.cost)]
    # update
    idsel = c(idsel, fopt)
    idcan = setdiff(idcan, fopt)
  }

  #  select and projection
  idxvec     = idsel
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = X%*%projection
  result$featidx = idxvec
  result$projection = projection
  result$algorithm  = "feature:procrustes"
  return(structure(result, class="Rdimtools"))
}



# auxiliary functions -----------------------------------------------------
#' @keywords internal
procrustes_cost <- function(Y, Z, k){
  if (is.vector(Y)){
    Y = matrix(Y, ncol=1)
  }
  if (is.vector(Z)){
    Z = matrix(Z, ncol=1)
  }
  term1 = base::sum(base::diag(Y%*%t(Y)))
  term2 = base::sum(base::diag(Z%*%t(Z)))
  term3 = -2*base::sum(base::svd(t(Z)%*%Y)$d[1:k])
  return(term1+term2+term3)
}
