#' Nonnegative Orthogonal Locality Preserving Projection
#'
#'
#' @references
#' \insertRef{zafeiriou_nonnegative_2010}{Rdimtools}
#'
#' @rdname linear_NOLPP
#' @author Kisung You
#' @export
do.nolpp <- function(X, ndim=2, type=c("proportion",0.1),
                     preprocess=c("center","decorrelate","whiten"),
                     t=1.0, maxiter=1000, reltol=1e-5){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.nolpp : 'ndim' is a positive integer in [1,#(covariates)].")
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
  #   5. t = kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t, .Machine$double.eps, Inf, compact=TRUE)){stop("* do.nolpp : 't' should be a bandwidth parameter in (0,Inf).")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
  #   2. neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask
  #   3. Dsqmat with kernelization
  Dsqmat = exp(-(as.matrix(dist(pX))^2)/t)

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR NONNEGATIVE OLPP
  #   1. compute auxiliary matrices
  A = Dsqmat*nbdmask
  L = diag(rowSums(A))-A

  #   2. compute cost function
  C = t(pX)%*%L%*%pX

  #   3. set initial matrix
  Uinit = matrix(runif(p*ndim),nrow=p)

  #   4. solve the minimization problem
  projection = aux.adjprojection(method_nnprojmin(C, Uinit, reltol, maxiter))


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
