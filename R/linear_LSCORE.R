#' Laplacian Score
#'
#' Laplacian Score (LSCORE) is an unsupervised linear feature extraction method. For each
#' feature/variable, it computes Laplacian score based on an observation that data from the
#' same class are often close to each other. Its power of locality preserving property is used, and
#' the algorithm selects variables with smallest scores.
#'
#' @examples
#' ## use iris data
#' ## it is known that feature 3 and 4 are more important.
#' data(iris)
#' set.seed(100)
#' subid    <- sample(1:150, 50)
#' iris.dat <- as.matrix(iris[subid,1:4])
#' iris.lab <- as.factor(iris[subid,5])
#'
#' ## try different kernel bandwidth
#' out1 = do.lscore(iris.dat, t=0.1)
#' out2 = do.lscore(iris.dat, t=1)
#' out3 = do.lscore(iris.dat, t=10)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=iris.lab, main="bandwidth=0.1")
#' plot(out2$Y, pch=19, col=iris.lab, main="bandwidth=1")
#' plot(out3$Y, pch=19, col=iris.lab, main="bandwidth=10")
#' par(opar)
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
#' @param t bandwidth parameter for heat kernel in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{lscore}{a length-\eqn{p} vector of laplacian scores. Indices with smallest values are selected.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @references
#' \insertRef{he_laplacian_2005}{Rdimtools}
#'
#' @rdname linear_LSCORE
#' @author Kisung You
#' @concept linear_methods
#' @export
do.lscore <- function(X, ndim=2, type=c("proportion",0.1),
                      preprocess=c("null","center","scale","cscale","whiten","decorrelate"), t=10.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.lscore : 'ndim' is a positive integer in [1,#(covariates)].")
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
  #   5. t : kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t, 1e-15, Inf, compact=TRUE)){stop("* do.lscore : 't' is a kernel bandwidth parameter in (0,Inf).")}

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
  ## COMPUTATION : MAIN PART FOR LAPLACIAN SCORE
  #   1. weight matrix
  Dsqmat  = exp(-(as.matrix(dist(pX))^2)/t)
  S       = Dsqmat*nbdmask
  diag(S) = 0
  #   2. auxiliary matrices
  D = diag(rowSums(S))
  L = D-S
  #   3. compute Laplacian score
  n1 = as.vector(rep(1,n))
  D1 = as.vector(D%*%matrix(rep(1,n)))
  fscore = rep(0,p)
  for (j in 1:p){
    # 3-1. select each feature
    fr = as.vector(pX[,j])
    # 3-2. adjust fr
    corrector  = as.double(sum(fr*D1)/sum(n1*D1))
    frtilde    = fr-corrector
    matfrtilde = matrix(frtilde)
    # 3-3. compute the score
    term1 = sum(as.vector(L%*%matfrtilde)*frtilde)
    term2 = sum(as.vector(D%*%matfrtilde)*frtilde)
    fscore[j] = term1/term2
  }
  #   4. select the smallest ones
  idxvec = base::order(fscore, decreasing=FALSE)[1:ndim]
  #   5. find the projection matrix
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$lscore  = fscore
  result$featidx = idxvec
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
