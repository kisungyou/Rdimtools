#' Sample-Dependent Locality Preserving Projection
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param t kernel bandwidth in \eqn{(0,\infty)}.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center" and other options of "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#'
#'
#' @examples
#' ## generate data
#' X <-
#'
#'
#' @references
#' \insertRef{yang_sample-dependent_2010}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_SDLPP
#' @export
do.sdlpp <- function(X, ndim=2, t = 1.0,
                     preprocess=c("center","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.sdlpp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. t
  t = as.double(t)
  if (!check_NumMM(t,0,1e+10,compact=FALSE)){stop("* do.sdlpp : 't' should be a positive real number.")}
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. preprocessing of data matrix
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"

  #   2. preliminary computations
  #   2-1. compute D : normalized squared distance
  D = (as.matrix(dist(pX, method="euclidean"))^2)
  for (i in 1:n){
    sumvecD = sum(D[i,])
    D[i,] = D[i,]/sumvecD
  }
  #   2-2. compute Ss : exponentiated D squared                  :: Problem with 'zerodiag'.
  Ss = exp(-D/(2*(t^2)))
  # if (diagzero){                                               :: follow the method directly.
  #   diag(Ss) = 0
  # }
  #   2-3. compute Ws : conditionally
  Ws = array(0,c(n,n))
  rowMeansSs = rowMeans(Ss)
  for (i in 1:n){
    for (j in 1:n){
      if (Ss[i,j] > (rowMeansSs[i])){
        Ws[i,j] = Ss[i,j]
      }
    }
  }
  #   2-4. compute Ds and Ls
  Wtilde = Ws+t(Ws)
  Ds = diag(rowSums(Ws))+diag(colSums(Ws))
  Ls = Ds-Wtilde

  #   3. main LPP part
  LHS = t(pX)%*%Ls%*%pX
  RHS = t(pX)%*%Ds%*%pX

  #   4. compute Projection Matrix
  geigs = geigen::geigen(LHS, RHS, TRUE)
  projection = matrix(geigs$vectors[,1:ndim],nrow=p)
  eigenvalue = as.vector(geigs$values[1:ndim])

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$eigval = eigenvalue
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
