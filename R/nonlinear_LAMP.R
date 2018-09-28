#' Local Affine Multidimensional Projection
#'
#' Local Affine Mulditimensional Projection (\emph{LAMP}) can be considered as
#' a nonlinear method even though each datam is projected using locally estimated
#' affine mapping. It first finds a low-dimensional embedding for control points
#' and then locates the rest data using affine mapping. We use \eqn{\sqrt{n}} number
#' of data as controls and Stochastic Neighborhood Embedding is applied as an
#' initial projection of control set. Note that this belongs to the method for
#' visualization so projection onto \eqn{\mathbf{R}^2} is suggested for use.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \dontrun{
#' ## load iris data
#' data(iris)
#' X <- as.matrix(iris[,1:4])
#'
#' ## let's compare with PCA
#' out1 <- do.pca(X, ndim=2)      # PCA
#' out2 <- do.lamp(X, ndim=2)     # LAMP
#'
#' ## visualize
#' par(mfrow=c(1,2))
#' plot(out1$Y[,1], out1$Y[,2], main="PCA")
#' plot(out2$Y[,1], out2$Y[,2], main="LAMP")
#' }
#'
#' @references
#' \insertRef{joia_local_2011}{Rdimtools}
#'
#' @seealso \code{\link{do.sne}}
#' @author Kisung You
#' @rdname nonlinear_LAMP
#' @export
do.lamp <- function(X, ndim=2, preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  ########################################################################
  ## 1. Type Checking
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lamp : 'ndim' is a positive integer in [1,#(covariates)).")}
  algpreprocess = match.arg(preprocess)

  controls = sort(sample(1:n, round(sqrt(n)), replace=FALSE))
  ncontrol = length(controls)

  ########################################################################
  ## 2. Preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  pXcontrol = pX[controls,]
  pXrest    = pX[-controls,]

  ########################################################################
  ## 3. Preliminary Computation
  pYcontrol = suppressMessages(do.sne(pXcontrol, ndim=ndim, preprocess="null")$Y)
  pYrest    = array(0,c(n-ncontrol, ndim))

  ########################################################################
  ## 4. Main Computation
  for (i in 1:(n-ncontrol)){
    # 4-1. target x
    x = pXrest[i,]
    # 4-2. compute weights
    alpha = rep(0,ncontrol)
    for (j in 1:ncontrol){
      alpha[j] = 1/sum((as.vector(x)-as.vector(pXcontrol[j,]))^2)
    }
    alphasum = sum(alpha)
    # 4-3. xbar and ybar
    xbar = as.vector(rep(0,p))
    ybar = as.vector(rep(0,ndim))
    for (j in 1:ncontrol){
      xbar = xbar + alpha[j]*as.vector(pXcontrol[j,])
      ybar = ybar + alpha[j]*as.vector(pYcontrol[j,])
    }
    xbar = xbar/alphasum
    ybar = ybar/alphasum  # still they are vectors
    # 4-4. build A and B
    A = array(0,c(ncontrol,p))
    B = array(0,c(ncontrol,ndim))
    for (j in 1:ncontrol){
      A[j,] = sqrt(alpha[j])*(as.vector(pXcontrol[j,])-xbar)
      B[j,] = sqrt(alpha[j])*(as.vector(pYcontrol[j,])-ybar)
    }
    # 4-5. svd decomposition
    svdAB = svd((t(A)%*%B))
    M     = (svdAB$u %*% t(svdAB$v))
    # 4-6. record it
    pYrest[i,] = as.vector(matrix(x-xbar, nrow=1)%*%M) + ybar;
  }
  # 4-6. rearrange
  pY = array(0,c(n,ndim))
  pY[controls,]  = pYcontrol
  pY[-controls,] = pYrest


  ########################################################################
  ## 5. return output
  result = list()
  result$Y = pY
  trfinfo$algtype = "nonlinear"
  result$trfinfo  = trfinfo
  return(result)
}


