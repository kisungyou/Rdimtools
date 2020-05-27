#' Unsupervised Discriminative Features Selection
#'
#' Though it may sound weird, this method aims at finding discriminative features
#' under the unsupervised learning framework. It assumes that the class label
#' could be predicted by a linear classifier and iteratively updates its
#' discriminative nature while attaining row-sparsity scores for selecting features.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param lbd regularization parameter for local Gram matrix to be invertible.
#' @param gamma regularization parameter for row-sparsity via \eqn{\ell_{2,1}} norm.
#' @param k size of nearest neighborhood for each data point.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
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
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150, 50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' #### try different neighborhood size
#' out1 = do.udfs(X, k=5)
#' out2 = do.udfs(X, k=10)
#' out3 = do.udfs(X, k=25)
#'
#' #### visualize
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=label, main="UDFS::k=5")
#' plot(out2$Y, pch=19, col=label, main="UDFS::k=10")
#' plot(out3$Y, pch=19, col=label, main="UDFS::k=25")
#' par(opar)
#'
#' @references
#' \insertRef{yang_l2_2011}{Rdimtools}
#'
#' @author Kisung You
#' @rdname feature_UDFS
#' @concept feature_methods
#' @export
do.udfs <- function(X, ndim=2, lbd=1.0, gamma=1.0, k=5,
                    preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.udfs : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. lbd and gamma
  if ((length(lbd)>1)||(lbd<0)){
    stop("* do.udfs : 'lbd' should be a nonnegative real number.")
  }
  if ((length(gamma)>1)||(gamma<0)){
    stop("* do.udfs : 'gamma' should be a nonnegative real number.")
  }
  #   5. k : neighborhood size
  k = round(k)
  if ((length(k)>1)||(k<1)||(k>=(nrow(X)-1))){
    stop("* do.udfs : 'k' should be a positive integer in [1,nrow(X)-2].")
  }

  #------------------------------------------------------------------------
  ## COMPUTATION
  #   1. data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  D       = as.matrix(dist(pX))

  #   2. sub-R routine : compute M
  M       = do.udfs.pre(t(pX), D, lbd, k)
  #   3. sub-R routine : score
  score   = do.udfs.score(M, ndim, gamma)
  idxvec  = base::order(score, decreasing=TRUE)[1:ndim]
  #   4. get the projection based on top scored ones
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

#' @keywords internal
#' @noRd
do.udfs.pre <- function(X, D, lambda, k){ # column stacking so be careful
  d = nrow(X) # column stacking
  n = ncol(X)

  Hk1 = diag(rep(1,k+1)) - (1/(k+1))*array(1,c(k+1,k+1))
  Dk1 = diag(rep(1,k+1))
  M   = array(0,c(n,n))
  for (i in 1:n){
    ordi = order(D[i,])[1:(k+1)] # from the front, use (k+1)
    Xi   = X[,ordi]
    Si   = array(0,c(n,(k+1)))
    for (p in 1:n){
      for (q in 1:(k+1)){
        if (p==ordi[q]){
          Si[p,q] = 1
        }
      }
    }
    XiS = (Xi%*%Hk1)
    Bi  = aux.pinv(t(XiS)%*%(XiS) + (lambda*Dk1))
    Mi  = (Si%*%Hk1%*%Bi%*%Hk1%*%t(Si))
    M   = (M + Mi)
  }
  outM = (X%*%M%*%t(X))
}

#' @keywords internal
#' @noRd
do.udfs.score <- function(M, ndim, gamma){
  d    = nrow(M)
  Dold = diag(rep(1,d))
  Dnew = array(0,c(d,d))

  increment = 100000
  maxiter   = 496
  citer     = 1
  verysm    = sqrt(10*.Machine$double.eps)
  while (increment > 1e-6){
    P = M + gamma*Dold
    W = RSpectra::eigs_sym(P, ndim, which="SM")$vectors

    for (i in 1:d){
      Dnew[i,i] = 1/max(c(2*sqrt(sum(as.vector(W[i,])^2)), verysm))
    }
    increment = base::norm(Dold-Dnew,"F")

    ## update
    Dold  = Dnew
    citer = citer + 1
    if (citer>=maxiter){
      break
    }
  }
  score = rep(0,d)
  for (i in 1:d){
    score[i] = sqrt(sum(as.vector(W[i,])^2))
  }
  return(score)
}
