#' Local Linear Laplacian Eigenmaps
#'
#' Local Linear Laplacian Eigenmaps is an unsupervised manifold learning method as an
#' extension of Local Linear Embedding (\code{\link{do.lle}}). It is claimed to be
#' more robust to local structure and noises. It involves the concept of
#' artificial neighborhood in constructing the adjacency graph for reconstruction of
#' the approximated manifold.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is \code{"null"}. See also \code{\link{aux.preprocess}} for more details.
#' @param K size of near neighborhood for each data point.
#' @param P size of artifical neighborhood.
#' @param bandwidth scale parameter for Gaussian kernel. It should be in \eqn{(0,1)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \dontrun{
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.integer(iris$Species)
#'
#' # see the effect bandwidth
#' out1 = do.llle(X, bandwidth=0.1)
#' out2 = do.llle(X, bandwidth=0.5)
#' out3 = do.llle(X, bandwidth=0.9)
#'
#' # visualize the results
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(out1$Y, col=label, main="bandwidth=0.1")
#' plot(out2$Y, col=label, main="bandwidth=0.5")
#' plot(out3$Y, col=label, main="bandwidth=0.9")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{liu_local_2016}{Rdimtools}
#'
#' @rdname nonlinear_LLLE
#' @seealso \code{\link{do.lle}}
#' @author Kisung You
#' @export
do.llle <- function(X, ndim=2, preprocess=c("null","center","scale","cscale","decorrelate","whiten"),
                    K=round(nrow(X)/2), P=max(round(nrow(X)/4),2), bandwidth=0.2){
  #------------------------------------------------------------------------
  ## PREPROCESSING : PARAMETERS
  #   1. data matrix
  aux.typecheck(X)
  N = nrow(X)
  #   2. ndim and neigs
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.llle : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  #   3. preprocessing of data
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. K the size of neighbors
  if ((length(K)>1)||(K<=1)||(K>=(N-1))){
    stop("* do.llle : 'K' must be in [2,N-1).")
  }
  #   5. P the size of artificial nbd
  if ((length(P)>1)||(P<=1)||(P>=K)){
    stop("* do.llle : 'P' must be in [2,K).")
  }
  #   6. bandwidth
  if ((length(bandwidth)>1)||(bandwidth<=0)||(bandwidth>=1)){
    stop("* do.llle : 'bandwidth' must be in (0,1).")
  }
  K = round(K)
  P = round(P)


  #------------------------------------------------------------------------
  ## PREPROCESSING : DATA
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  D       = as.matrix(stats::dist(pX))
  for (i in 1:N){
    D[i,i] = Inf
  }

  #------------------------------------------------------------------------
  ## STEP 1 : NBD IDENTIFICATION
  di.vec = rep(0,N)
  Ci.vec = list()
  Ii.vec = list()
  for (i in 1:N){
    tgtD = as.vector(D[i,])
    idxK = which(tgtD <= sort(tgtD)[K])
    if (length(idxK)>K){
      idxK = sample(idxK, K, replace=FALSE)
    }


    Zi = pX[idxK,]-matrix(rep(as.vector(pX[i,]),K),nrow=K,byrow=TRUE)
    Ci.vec[[i]] = Zi%*%t(Zi)   # cmopute K-by-K matrix
    Ii.vec[[i]] = idxK         # index of K neighbors
    di.vec[i]   = sum(diag(Ci.vec[[i]]))
  }

  #------------------------------------------------------------------------
  ## STEP 2 : MAIN COMPUTATION
  #   2-1. prepare empty arrays
  out.D = array(0,c(N,N))
  out.L = array(0,c(N,N))
  #   2-2. iterate over each point
  for (i in 1:N){
    # 1. prepare Ai which will be used over, and others
    if (P>1){
      Ai.mat = rbind(rep(1,K),array(0,c(P-1,K))) # Ai.mat = PxK
    } else {
      Ai.mat = matrix(rep(1,K),nrow=1)
    }
    Ci     = Ci.vec[[i]]
    Ci.inv = aux.pinv(Ci)
    di     = di.vec[i]
    lbdi   = rep(0,P)
    # 2. iterate over p
    for (p in 1:P){
      if (p==1){
        Ai = Ai.mat[1,]
      } else {
        Ai = (t(Ai.mat[1:p,]))
      }
      wi = as.vector(Ci.inv%*%Ai%*%solve(t(Ai)%*%Ci.inv%*%Ai, llle.bp(p))) # should be (Kx1) but it's vector now
      if (p<P){
        Ai.mat[p+1,] = wi
      }
      lbdi[p] = sum(as.vector(Ci%*%matrix(wi,ncol=1))*wi)
    }
    # 3. Hi
    Hi = array(0,c(K,P))
    for (p in 1:P){
      Hi[,p] = exp(-lbdi[p]/(2*bandwidth*di))*as.vector(Ai.mat[p,])
    }
    si = as.vector(colSums(Hi))
    # 4. update
    Ii = Ii.vec[[i]]
    out.D[i,i]   = sum(si*si)
    out.L[i,i]   = out.L[i,i]   + sum(si*si)
    out.L[Ii,i]  = out.L[Ii,i]  - Hi%*%matrix(si,ncol=1)
    out.L[i,Ii]  = out.L[i,Ii]  - matrix(si,nrow=1)%*%t(Hi)
    out.L[Ii,Ii] = out.L[Ii,Ii] + Hi%*%t(Hi)
  }
  #   2-3. geigen
  Ylow = aux.geigen(out.L, out.D, (ndim+1), maximal=FALSE)

  #------------------------------------------------------------------------
  ## RETURN
  result   = list()
  result$Y = Ylow[,2:(ndim+1)]
  result$trfinfo = trfinfo
  return(result)
}

#' @keywords internal
#' @noRd
llle.bp <- function(p){
  if (p==1){
    return(matrix(1,nrow=1))
  } else {
    return(rbind(1, array(0,c(p-1,1))))
  }
}
