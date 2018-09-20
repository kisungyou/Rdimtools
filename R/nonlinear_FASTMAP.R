#' FastMap
#'
#' \code{do.fastmap} is an implementation of \emph{FastMap} algorithm. Though
#' it shares similarities with MDS, it is innately a nonlinear method that iteratively updates
#' the projection information using pairwise distance information.
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
#' ## let's compare with other methods
#' out1 <- do.pca(X, ndim=2)      # PCA
#' out2 <- do.mds(X, ndim=2)      # Classical MDS
#' out3 <- do.fastmap(X, ndim=2)  # FastMap
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="PCA")
#' plot(out2$Y[,1], out2$Y[,2], main="MDS")
#' plot(out3$Y[,1], out3$Y[,2], main="FastMap")
#' }
#'
#' @references
#' \insertRef{faloutsos_fastmap:_1995}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_FastMap
#' @export
do.fastmap <- function(X, ndim=2, preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  ########################################################################
  ## 1. Type Checking
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.fastmap : 'ndim' is a positive integer in [1,#(covariates)).")}
  algpreprocess = match.arg(preprocess)
  k = ndim

  ########################################################################
  ## 2. Preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  ########################################################################
  ## 3. Main Computation
  Dold    = as.matrix(dist(pX))
  Dnew    = array(0,c(n,n))
  output  = array(0,c(n,ndim))
  for (i in 1:k){
    # 3-1. find row and column index for maximal element
    maxidx = aux.findmaxidx(Dold)
    ida    = maxidx[1]
    idb    = maxidx[2]
    # 3-2. precompute some values
    dab2  = (sum(as.vector(pX[ida,]-pX[idb,])^2))
    if (dab2 > 10*.Machine$double.eps){
      dab   = sqrt(dab2)
      # 3-3. compute coefficient
      for (j in 1:n){
        output[j,i] = (sum(as.vector(pX[ida,]-pX[j,])^2) + dab2 - sum(as.vector(pX[idb,]-pX[j,])^2))/(2*dab)
      }
    } # or, leave it as zero

    # 3-4. update D : compute and alter
    for (it1 in 1:n){
      for (it2 in 1:n){
        Dnew[it1,it2] = sqrt((Dold[it1,it2]^2) - ((output[it1,i]-output[it2,i])^2))
      }
    }
    Dold = Dnew

  }

  ########################################################################
  ## 4. return output
  result = list()
  result$Y = output
  trfinfo$algtype = "nonlinear"
  result$trfinfo  = trfinfo
  return(result)
}
