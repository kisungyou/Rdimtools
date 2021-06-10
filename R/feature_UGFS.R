#' Unsupervised Graph-based Feature Selection
#'
#' UGFS is an unsupervised feature selection method with two parameters \code{nbdk} and \code{varthr} that it constructs
#' an affinity graph using local variance computation and scores variables based on PageRank algorithm.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param nbdk the size of neighborhood for local variance computation.
#' @param varthr threshold value for affinity graph construction. If too small so that the graph of variables is not constructed, it returns an error.
#' @param preprocess an additional option for preprocessing the data. Default is "null". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{prscore}{a length-\eqn{p} vector of score computed from PageRank algorithm. Indices with largest values are selected.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' ## it is known that feature 3 and 4 are more important.
#' data(iris)
#' iris.dat <- as.matrix(iris[,1:4])
#' iris.lab <- as.factor(iris[,5])
#'
#' ## try multiple thresholding values
#' out1 = do.ugfs(iris.dat, nbdk=10, varthr=0.5)
#' out2 = do.ugfs(iris.dat, nbdk=10, varthr=5.0)
#' out3 = do.ugfs(iris.dat, nbdk=10, varthr=9.5)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=iris.lab, main="bandwidth=0.1")
#' plot(out2$Y, pch=19, col=iris.lab, main="bandwidth=1")
#' plot(out3$Y, pch=19, col=iris.lab, main="bandwidth=10")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{henni_unsupervised_2018}{Rdimtools}
#'
#' @rdname feature_UGFS
#' @author Kisung You
#' @concept feature_methods
#' @export
do.ugfs <- function(X, ndim=2, nbdk=5, varthr=2.0, preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.ugfs : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. extra parameters
  myk   = round(nbdk)
  mythr = as.double(varthr)

  #------------------------------------------------------------------------
  ## COMPUTATION : Preliminary
  #  Preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #  Nearest-neighbor
  knnnbd  = RANN::nn2(pX, k=(myk+1))

  #------------------------------------------------------------------------
  ## COMPUTATION : Main
  #  Step 1. compute local variance
  Varji = array(0,c(n,p))
  for (i in 1:n){
    tgtvec = pX[i,]
    tgtmat = pX[as.vector(knnnbd$nn.idx[,2:(myk+1)]),]
    Varji[i,] = ugfs.var(tgtvec, tgtmat)
  }
  print(paste0("* do.ugfs : range of variances is [",round(min(Varji),2),",",round(max(Varji),2),"]."))

  #  Step 2. binarize
  VarBin = array(0,c(n,p))
  VarBin[(Varji <= mythr)] = 1
  if (all(VarBin==0)){
    stop("* do.ugfs : 'varthr' is too small that graph may not be constructed. Use larger value.")
  }
  #  Step 3. graph construction
  A = array(0,c(p,p))
  for (i in 1:n){
    Sp = c()
    for (j in 1:p){
      if (VarBin[i,j]>0.5){
        Sp = c(Sp, j)
      }
    }
    A[Sp,Sp] = 1
  }
  diag(A) = 0
  #  Step 4. compute the score via pagerank algorithm
  pgscore = as.vector(v2aux_pagerank(A))
  # pgscore = mypagerank(A)
  #  Step 5. select the largest ones and find the projection
  idxvec     = base::order(pgscore, decreasing=TRUE)[1:ndim]
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$A = A
  result$prscore = pgscore
  result$featidx = idxvec
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}



# auxiliary functions -----------------------------------------------------
#' @keywords internal
ugfs.var <- function(vec, mat){ # compute for all variables at once
  if (is.vector(mat)){
    mat = matrix(mat, nrow=1)
  }
  k = base::nrow(mat)
  p = base::ncol(mat)
  output = rep(0,p)
  for (i in 1:p){
    tgtval = vec[i]
    tgtvec = mat[,i]
    output[i] = sum((tgtval-tgtvec)^2)/k
  }
  return(output)
}
#' @keywords internal
mypagerank <- function(A){
  N = nrow(A)
  k = base::rowSums(A)
  M = t(diag(1/k)%*%A)
  d = 0.85

  Rold = rep(1/N,N)
  for (i in 1:100){
    Rnew = d*as.vector(M%*%Rold) + ((1-d)/N)*rep(1,N)
    Rinc = sqrt(sum((Rold-Rnew)^2))
    Rold = Rnew
    if (Rinc < 1e-6){
      break
    }
  }
  return(Rold)
}
