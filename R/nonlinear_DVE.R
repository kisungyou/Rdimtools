#' Distinguishing Variance Embedding
#'
#' Distinguishing Variance Embedding (DVE) is an unsupervised nonlinear manifold learning method.
#' It can be considered as a balancing method between Maximum Variance Unfolding and Laplacian
#' Eigenmaps. The algorithm unfolds the data by maximizing the global variance subject to the
#' locality-preserving constraint. Instead of defining certain kernel, it applies local scaling scheme
#' in that it automatically computes adaptive neighborhood-based kernel bandwidth.
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
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \donttest{
#' ## generate swiss-roll dataset of size 100
#' set.seed(100)
#' X <- aux.gensamples(dname="crown", n=100)
#'
#' ## try different nbd size
#' out1 <- do.dve(X, type=c("proportion",0.5))
#' out2 <- do.dve(X, type=c("proportion",0.7))
#' out3 <- do.dve(X, type=c("proportion",0.9))
#' }
#'
#' @references
#' \insertRef{wang_combining_2009}{Rdimtools}
#'
#' \insertRef{qinggang_distinguishing_2010}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_DVE
#' @concept nonlinear_methods
#' @export
do.dve <- function(X, ndim=2, type=c("proportion",0.1),
                   preprocess=c("null","center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.dve : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. type
  nbdtype = type
  nbdsymmetric = "union"
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : DATA PREPROCESSING
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN STOPS FOR DVE
  #   1. construct neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  G      = nbdmask
  Gprime = (!nbdmask)
  diag(Gprime)=FALSE   # adjust diagonal elements

    #   2. local scaling method
  locstd = rep(0,n)
  for (i in 1:n){
    #   2-1. selection of target vector and target matrix
    idxnbd = which(nbdmask[i,])
    tgtvec = as.vector(pX[i,])
    tgtmat = pX[idxnbd,]

    #   2-2. compute
    locstd[i] = dve_localscaling(tgtvec, tgtmat)
  }

  #   3. weights
  #   3-1. compute pairwise Euclidean distance and adjust it
  Wsqmat =  exp(-diag(1/locstd)%*%(as.matrix(dist(pX))^2)%*%diag(1/locstd))
  #   3-2. separate out
  W      = Wsqmat*G
  Wprime = Gprime*1.0
  #   3-3. compute auxiliary variables
  L      = diag(rowSums(W))-W
  Lprime = diag(rowSums(Wprime))-Wprime

  #   4. solve for geigen problem : use largest ones
  Youtput = aux.geigen(Lprime, L, ndim, maximal=TRUE)


  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN STOPS FOR DVE
  result = list()
  result$Y = Youtput
  result$trfinfo  = trfinfo
  return(result)
}





# . -----------------------------------------------------------------------
#' @keywords internal
#' @noRd
dve_localscaling <- function(vec, mat){
  if (is.vector(mat)){
    vecdiff = as.vector(vec)-as.vector(mat)
    result = sqrt(sum(vecdiff*vecdiff))
  } else {
    n   = nrow(mat)
    record = rep(0,n)
    for (i in 1:n){
      vecdiff = vec - as.vector(mat[i,])
      record[i] = sum(vecdiff*vecdiff)
    }
    result = sqrt(max(record))
    if (result<.Machine$double.eps){
      result = .Machine$double.eps
    }
  }
  return(result)
}
