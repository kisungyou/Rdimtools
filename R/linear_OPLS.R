#' Orthogonal Partial Least Squares
#'
#' Also known as multilinear regression or semipenalized CCA, Orthogonal Partial Least Squares (OPLS)
#' was first used to perform multilinear ordinary least squares. In its usage, unlike PLS or CCA,
#' OPLS does not rely on projected variance of response -or, \code{data2}. Instead, it exploits projected
#' variance of input - covariance of \code{data1} and relates it under cross-covariance setting. Therefore,
#' OPLS only returns projection information of \code{data1}, just like any other unsupervised methods in our package.
#'
#' @param data1 an \eqn{(n\times N)} data matrix whose rows are observations.
#' @param data2 an \eqn{(n\times M)} data matrix whose rows are observations.
#' @param ndim an integer-valued target dimension.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix of projected observations from \code{data1}.}
#' \item{projection}{an \eqn{(N\times ndim)} whose columns are loadings for \code{data1}.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction for \code{data1}.}
#' \item{eigvals}{a vector of eigenvalues for iterative decomposition.}
#' }
#'
#' @examples
#' ## generate 2 normal data matrices
#' mat1 = matrix(rnorm(100*12),nrow=100)+10 # 12-dim normal
#' mat2 = matrix(rnorm(100*6), nrow=100)-10 # 6-dim normal
#'
#' ## compare OPLS and PLS
#' res_opls = do.opls(mat1, mat2, ndim=2)
#' res_pls  = do.pls(mat1, mat2, ndim=2)
#'
#' ## visualize
#' par(mfrow=c(1,2))
#' plot(res_opls$Y[,1], res_opls$Y[,2], main="OPLS result")
#' plot(res_pls$Y1[,1], res_pls$Y1[,2], main="PLS result")
#'
#' @references
#' \insertRef{barker_partial_2003}{Rdimtools}
#'
#' @seealso \code{\link{do.pls}}
#' @author Kisung You
#' @rdname linear_OPLS
#' @export
do.opls <- function(data1,data2,ndim=2){
  ## Preprocessing
  #   1. data matrices type
  if (is.vector(data1)&&is.vector(data2)){stop("* do.opls : either one of 'data1' and 'data2' should be a matrix.")}
  if (is.vector(data1)){data1 = as.matrix(data1)}
  if (is.vector(data2)){data2 = as.matrix(data2)}
  aux.typecheck(data1)
  aux.typecheck(data2)
  #   2. size arguments
  if (nrow(data1)!=nrow(data2)){stop("* do.opls : both data matrices should have same number of observations.")}
  n = nrow(data1)
  N = ncol(data1)
  M = ncol(data2)
  #   3. dimension arguments
  if ((!check_ndim(ndim,N))||(!check_ndim(ndim,M))){
    stop("* do.opls : 'ndim' is a positive integer in [1,min(ncol(X),ncol(Y)).")
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
  xP = array(0,c(N,ndim)) # be careful about transpose

    #   2. iterations
  evals = rep(0,ndim)
  for (i in 1:ndim){
    Cx = (t(pX)%*%pX)/n
    Cy = (t(pY)%*%pY)/n
    Cxy= (t(pX)%*%pY)/n

    #   2-1. eigen decomposition
    LHS = Cxy%*%t(Cxy)
    RHS = Cx
    SOL = aux.bicgstab(RHS, LHS, verbose=FALSE)
    res = RSpectra::eigs_sym(SOL$x, 1, which="LA")
    evals[i] = res$values[1]

    tmpP = matrix(res$vectors)

    #   2-2. record projections and loadings
    xT[,i] = pX%*%tmpP
    xP[,i] = tmpP

    #   2-3. regress out!
    pX = aux_regout(pX, as.vector(tmpP))
  }

  ## Return output
  result = list()
  result$Y = xT
  result$projection = aux.adjprojection(xP)
  result$trfinfo    = trfinfoX
  result$eigvals       = evals
  return(result)
}
