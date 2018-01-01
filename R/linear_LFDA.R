#' Local Fisher Discriminant Analysis
#'
#'
#'
#' @references
#' \insertRef{sugiyama_local_2006}{Rdimtools}
#'
#' \insertRef{zelnik-manor_self-tuning_2005}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LFDA
#' @export
do.lfda <- function(X, label, ndim=2, preprocess=c("center","decorrelate","whiten"),
                    type=c("proportion",0.1), symmetric=c("union","intersect","asymmetric"),
                    localscaling=TRUE){
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
  for (i in 1:length(ulabel)){
    if (sum(label==ulabel[i])==1){
      stop("* do.lfda : no degerate class of size 1 is allowed.")
    }
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    warning("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lfda : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. nbd-type
  nbdtype = type
  #   6. nbd-symmetric
  if (missing(symmetric)){
    nbdsymmetric = "union"
  } else {
    nbdsymmetric = match.arg(symmetric)
  }
  #   7. localscaling
  if (!is.logical(localscaling)){
    stop("* do.lfda : 'localscaling' must be a logical flag.")
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. Preprocessing the data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
  #   2. neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  D     = nbdstruct$dist
  Dmask = nbdstruct$mask

}
