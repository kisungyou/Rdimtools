#' Nonnegative Orthogonal Neighborhood Preserving Projections
#'
#' @references
#' \insertRef{zafeiriou_nonnegative_2010}{Rdimtools}
#'
#' @rdname linear_NONPP
#' @author Kisung You
#' @export
do.nonpp <- function(X, ndim=2, type=c("proportion",0.1),
                     preprocess=c("center","decorrelate","whiten"),
                     maxiter=1000, reltol=1e-5){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.nonpp : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   3. type
  nbdtype = type
  nbdsymmetric = "union"
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY and LLE step
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"

  #   2. neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)

  #   3. LLE computation
  regparam = 1.0
  W = array(0,c(n,n))
  for (i in 1:n){
    #   3-1. separate target mask vector
    tgtidx  = which(nbdstruct$mask[i,])
    #   3-2. select data
    #        For convenience, target matrix is transposed for Armadillo
    vec_tgt = pX[i,]
    mat_tgt = t(pX[tgtidx,])
    k = ncol(mat_tgt)
    #   3-3. no automatic regularization
    W[i,tgtidx] = method_lleW(mat_tgt,vec_tgt,regparam);
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION FOR N-ONPP
  #   1. cost function : C
  diagN  = diag(n)
  C      = t(pX)%*%(t(diagN-W)%*%(diagN-W))%*%pX
  #   2. initialize for U
  Uinit = matrix(runif(p*ndim),nrow=p)
  #   3. compute projection matrix
  projection = aux.adjprojection(method_nnprojmin(C, Uinit, reltol, maxiter))

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
