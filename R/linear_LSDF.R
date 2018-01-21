#' Locality Sensitive Discriminant Feature
#'
#'
#' @rdname linear_LSDF
#' @author Kisung You
#' @export
do.lsdf <- function(X, label, ndim=2, type=c("proportion",0.1),
                    preprocess=c("null","center","whiten","decorrelate"), gamma=100){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label : check and return a de-factored vector
  #   For this example, there should be no degenerate class of size 1.
  label  = check_label(label, n)
  ulabel = unique(label)
  if (all(!is.na(ulabel))){
    message("* Semi-Supervised Learning : there is no missing labels. Consider using Supervised methods.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.lsdf : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. type
  nbdtype = type
  nbdsymmetric = "union"
  #   5. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   6. gamma
  gamma = as.double(gamma)
  if (!check_NumMM(gamma,1,1e+10)){stop("* do.lsdf : 'gamma' is a large positive real number.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  if (algpreprocess=="null"){
    trfinfo = list()
    trfinfo$type = "null"
    pX = X
  } else {
    tmplist = aux.preprocess(X,type=algpreprocess)
    trfinfo = tmplist$info
    pX      = tmplist$pX
  }
  trfinfo$algtype = "linear"

  #   2. build neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION FOR LSDF
  #   1. build Within-Class weight
  Sw = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
    }
  }
}
