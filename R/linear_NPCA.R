#' Nonnegative Principal Component Analysis
#'
#' @references
#' \insertRef{zafeiriou_nonnegative_2010}{Rdimtools}
#'
#' @rdname linear_NPCA
#' @author Kisung You
#' @export
do.npca <- function(X, ndim=2, preprocess=c("center","decorrelate","whiten"),
                    maxiter=1000, reltol=1e-5){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.npca : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : DATA PREPROCESSING
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR NONNEGATIVE PCA
  #   1. initialize for U
  Uinit = matrix(runif(p*ndim),nrow=p)
  #   2. compute C
  C = cov(pX)
  #   3. compute projection matrix
  projection = aux.adjprojection(method_nnprojmax(C, Uinit, reltol, maxiter))


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
