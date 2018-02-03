#' Fisher Score
#'
#' Fisher Score (FSCORE) is a supervised linear feature extraction method. For each
#' feature/variable, it computes Fisher score, a ratio of between-class variance to within-class variance.
#' The algorithm selects variables with largest Fisher scores and returns an indicator projection matrix.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null" and other options of "center", "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
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
#' \dontrun{
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=33)-100
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+100
#'
#' ## merge the data and create a label correspondingly
#' X      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## compare Fisher score with LDA
#' out1 = do.lda(X, label)
#' out2 = do.fscore(X, label)
#'
#' ## visualize
#' par(mfrow=c(1,2))
#' plot(out1$Y[,1], out1$Y[,2], main="LDA")
#' plot(out2$Y[,1], out2$Y[,2], main="Fisher Score")
#' }
#'
#' @references
#' \insertRef{fisher_use_1936}{Rdimtools}
#'
#' @rdname linear_FSCORE
#' @author Kisung You
#' @export
do.fscore <- function(X, label, ndim=2, preprocess=c("null","center","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label vector
  label  = check_label(label, n)
  ulabel = unique(label)
  C      = length(ulabel)
  if (C==1){
    stop("* do.fscore : 'label' should have at least 2 unique labelings.")
  }
  if (C==n){
    stop("* do.fscore : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.fscore : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PREPROCESSING OF THE DATA
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

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION FOR LSDF
  #   1. compute Fisher score for each feature
  fscore = rep(0,p)
  for (i in 1:p){
    vecfr     = as.vector(pX[,i])
    fscore[i] = fscore_single(vecfr, label, ulabel, C)
  }
  #   2. select the largest ones
  idxvec = base::order(fscore, decreasing=TRUE)[1:ndim]
  #   3. find the projection matrix
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





#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
fscore_single <- function(vec, label, ulabel, C){
  if (length(vec)!=length(label)){
    stop("* fscore_single : vec and label.")
  }
  if (length(ulabel)!=C){
    stop("* fscore_single : ulabel and C.")
  }

  # 1. compute overall values
  all_mean = mean(vec)
  all_var  = var(vec)

  # 2. class-wise information
  term1 = 0.0 # numerator
  term2 = 0.0 # denominator
  for (i in 1:C){
    idxc = which(label==ulabel[i])
    vecc = vec[idxc]
    ni   = length(idxc)

    term1 = term1 + ni*((mean(vecc)-all_mean)^2)
    term2 = term2 + ni*(var(vecc))
  }
  output = term1/term2
  return(output)
}
