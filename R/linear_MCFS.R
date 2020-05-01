#' Multi-Cluster Feature Selection
#'
#' Multi-Cluster Feature Selection (MCFS) is an unsupervised feature selection method. Based on
#' a multi-cluster assumption, it aims at finding meaningful features using sparse reconstruction of
#' spectral basis using LASSO.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param K assumed number of clusters in the original dataset.
#' @param lambda \eqn{\ell_1} regularization parameter in \eqn{(0,\infty)}.
#' @param t bandwidth parameter for heat kernel in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=20)-100
#' dt2  = aux.gensamples(n=20)
#' dt3  = aux.gensamples(n=20)+100
#'
#' ## merge the data and create a label correspondingly
#' X      = rbind(dt1,dt2,dt3)
#' label  = rep(1:3, each=20)
#'
#' ## try different regularization parameters
#' out1 = do.mcfs(X, lambda=0.01)
#' out2 = do.mcfs(X, lambda=0.1)
#' out3 = do.mcfs(X, lambda=1)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=label, main="lambda=0.01")
#' plot(out2$Y, pch=19, col=label, main="lambda=0.1")
#' plot(out3$Y, pch=19, col=label, main="lambda=1")
#' par(opar)
#'
#' @references
#' \insertRef{cai_unsupervised_2010}{Rdimtools}
#'
#' @rdname linear_MCFS
#' @author Kisung You
#' @concept linear_methods
#' @export
do.mcfs <- function(X, ndim=2, type=c("proportion",0.1),
                    preprocess=c("null","center","scale","cscale","whiten","decorrelate"),
                    K=max(round(nrow(X)/5),2), lambda=1.0, t=10.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.mcfs : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   3. type
  nbdtype = type
  nbdsymmetric = "union"
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. K : cluster numbers
  K = as.integer(K)
  if (!check_NumMM(K,2,round(nrow(X)/5))){stop("* do.mcfs : 'K' is an assumed cluster size in [2,#(samples)/2].")}
  #   6. lambda
  lambdaval = as.double(lambda)
  if (!check_NumMM(lambdaval,0,Inf,compact=FALSE)){stop("* do.mcfs : 'lambda' is a LASSO parameter in (0,Inf).")}
  #   7. t : bandwidth parameter
  t = as.double(t)
  if (!check_NumMM(t,1e-10,Inf,compact=TRUE)){stop("* do.mcfs : 't' is a bandwidth parameter in (0,Inf).")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. build neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR MULTI-CLUSTER FEATURE SELECTION
  #   1. construct nbd graph with weights; W
  Dsqmat = exp(-(as.matrix(dist(pX))^2)/t)
  W      = Dsqmat*nbdmask
  #   2. solve generalized eigenvalue problem
  D = diag(rowSums(W))
  L = D-W
  Y = aux.geigen(L,D,K,maximal=FALSE)
  #   3. solve K number of LASSO problems
  A = array(0,c(p,K))
  for (i in 1:K){
    # 3-1. take one column vector
    y      = as.vector(Y[,i])
    # 3-2. solve with LASSO; I will do it with mine
    # solved = ADMM ::admm.lasso(pX, y, lambda=lambdaval)$x
    # 3-3. record the solved
    # A[,i] = as.vector(solved)
    A[,i] = as.vector(admm_lasso(pX, y, lambdaval))
  }
  #   4. find the solution
  fscore = rep(0,p)
  for (i in 1:p){
    # 4-1. select one vector
    veca = base::abs(as.vector(A[i,]))
    # 4-2. take the largest
    fscore[i] = max(veca)
  }
  #   5. select the largest ones
  idxvec = base::order(fscore, decreasing=TRUE)[1:ndim]
  #   6. find the projection matrix
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$featidx = idxvec
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
