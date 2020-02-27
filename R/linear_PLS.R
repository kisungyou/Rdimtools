#' Partial Least Squares
#'
#' Given two data sets, Partial Least Squares (PLS) aims at maximizing cross-covariance of latent variables for each data matrix,
#' therefore it can be considered as supervised methods. As we have two input matrices, \code{do.pls} generates two sets of
#' outputs. Though it is widely used for regression problem, we used it in dimension reduction setting. For
#' algorithm aspects, we used recursive gram-schmidt orthogonalization in conjunction with extracting projection vectors under
#' eigen-decomposition formulation, as the problem dimension matters only up to original dimensionality.
#' For more details, see \href{https://en.wikipedia.org/wiki/Partial_least_squares_regression}{Wikipedia entry} on PLS.
#'
#' @param data1 an \eqn{(n\times N)} data matrix whose rows are observations
#' @param data2 an \eqn{(n\times M)} data matrix whose rows are observations
#' @param ndim an integer-valued target dimension.
#'
#' @return a named list containing
#' \describe{
#' \item{Y1}{an \eqn{(n\times ndim)} matrix of projected observations from \code{data1}.}
#' \item{Y2}{an \eqn{(n\times ndim)} matrix of projected observations from \code{data2}.}
#' \item{projection1}{an \eqn{(N\times ndim)} whose columns are loadings for \code{data1}.}
#' \item{projection2}{an \eqn{(M\times ndim)} whose columns are loadings for \code{data2}.}
#' \item{trfinfo1}{a list containing information for out-of-sample prediction for \code{data1}.}
#' \item{trfinfo2}{a list containing information for out-of-sample prediction for \code{data2}.}
#' \item{eigvals}{a vector of eigenvalues for iterative decomposition.}
#' }
#'
#' @examples
#' ## generate 2 normal data matrices
#' mat1 = matrix(rnorm(100*12),nrow=100)+10 # 12-dim normal
#' mat2 = matrix(rnorm(100*6), nrow=100)-10 # 6-dim normal
#'
#' ## project onto 2 dimensional space for each data
#' output = do.pls(mat1, mat2, ndim=2)
#'
#' ## visualize
#' opar <- par(mfrow=c(1,2), no.readonly=TRUE)
#' plot(output$Y1, main="proj(mat1)")
#' plot(output$Y2, main="proj(mat2)")
#' par(opar)
#'
#' @references
#' \insertRef{wold_path_1975}{Rdimtools}
#'
#' \insertRef{rosipal_overview_2006}{Rdimtools}
#'
#'
#' @seealso \code{\link{do.cca}}
#' @author Kisung You
#' @rdname linear_PLS
#' @export
do.pls <- function(data1,data2,ndim=2){
  ## Preprocessing
  #   1. data matrices type
  if (is.vector(data1)&&is.vector(data2)){stop("* do.pls : either one of 'data1' and 'data2' should be a matrix.")}
  if (is.vector(data1)){data1 = as.matrix(data1)}
  if (is.vector(data2)){data2 = as.matrix(data2)}
  aux.typecheck(data1)
  aux.typecheck(data2)
  #   2. size arguments
  if (nrow(data1)!=nrow(data2)){stop("* do.pls : both data matrices should have same number of observations.")}
  n = nrow(data1)
  N = ncol(data1)
  M = ncol(data2)
  #   3. dimension arguments
  if ((!check_ndim(ndim,N))||(!check_ndim(ndim,M))){
    stop("* do.pls : 'ndim' is a positive integer in [1,min(ncol(X),ncol(Y)).")
  }
  ndim = as.integer(ndim)
  #   4. perform centering
  tmplistX = aux.preprocess(data1,type="center")
  tmplistY = aux.preprocess(data2,type="center")

  trfinfoX = tmplistX$info; trfinfoX$algtype = "linear"
  trfinfoY = tmplistY$info; trfinfoY$algtype = "linear"

  pX = tmplistX$pX
  pY = tmplistY$pX

  ## Main Computation
  #   1. get ready for data generation process
  xT = array(0,c(n,ndim))
  yU = array(0,c(n,ndim))
  xP = array(0,c(N,ndim)) # be careful about transpose
  yQ = array(0,c(M,ndim)) # be careful about transpose

  zerosNN = array(0,c(N,N))
  zerosMM = array(0,c(M,M))
  #   2. iterations
  evals = rep(0,ndim)
  for (i in 1:ndim){
    Cx = (t(pX)%*%pX)/n
    Cy = (t(pY)%*%pY)/n
    Cxy= (t(pX)%*%pY)/n

    #   2-1. eigen decomposition
    LHS = rbind(cbind(zerosNN, Cxy), cbind(t(Cxy), zerosMM))
    res = RSpectra::eigs_sym(LHS, 1, which="LA")
    tmpP = matrix(res$vectors[1:N])
    tmpQ = matrix(res$vectors[(N+1):(N+M)])
    evals[i] = res$values[1]

    #   2-2. record projections and loadings
    xT[,i] = pX%*%tmpP
    yU[,i] = pY%*%tmpQ
    xP[,i] = tmpP
    yQ[,i] = tmpQ

    #   2-3. regress out!
    pX = aux_regout(pX, as.vector(tmpP))
    pY = aux_regout(pY, as.vector(tmpQ))
  }

  ## Return output
  result = list()
  result$Y1 = xT
  result$Y2 = yU
  result$projection1 = xP
  result$projection2 = yQ
  result$trfinfo1    = trfinfoX
  result$trfinfo2    = trfinfoY
  result$eigvals       = evals
  return(result)
}
