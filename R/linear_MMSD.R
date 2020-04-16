#' Multiple Maximum Scatter Difference
#'
#' Multiple Maximum Scatter Difference (MMSD) is a supervised linear dimension reduction method.
#' It is a variant of MSD in that discriminant vectors are orthonormal. Similar to MSD, it also does not suffer from
#' rank deficiency issue of scatter matrix.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param C nonnegative balancing parameter for intra- and inter-class scatter.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \donttest{
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=33)-100
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+100
#'
#' ## merge the data and create a label correspondingly
#' X      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## try different balancing parameter
#' out1 = do.mmsd(X, label, C=0.01)
#' out2 = do.mmsd(X, label, C=1)
#' out3 = do.mmsd(X, label, C=100)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="MMSD::C=0.01")
#' plot(out2$Y, main="MMSD::C=1")
#' plot(out3$Y, main="MMSD::C=100")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{fengxisong_multiple_2007}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_MMSD
#' @concept linear_methods 
#' @export
do.mmsd <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","whiten","decorrelate"), C=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label
  label  = as.numeric(as.factor(label))
  ulabel = unique(label)
  K      = length(ulabel)
  if (K==1){
    stop("* do.mmsd : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    stop("* do.mmsd : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.mmsd : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. C : nonnegative constant which balances two objective functions
  C = as.double(C)
  if (!check_NumMM(C, 0, 1e+10, compact=TRUE)){stop("* do.mmsd : 'C' is a nonnegative balancing parameter.")}



  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX


  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR MSD
  # 1. compute S_W  : within-group variance for multiclss problem
  Sw = array(0,c(p,p))
  for (i in 1:K){
    idxnow = which(label==ulabel[i])
    Sw     = Sw + lda_outer(pX[idxnow,])
  }
  # 2. compute S_B  : between-group variance for multiclass problem
  Sb = array(0,c(p,p))
  m  = colMeans(pX)
  for (i in 1:K){
    idxnow = which(label==ulabel[i])
    Nk     = length(idxnow)
    mdiff  = (colMeans(pX[idxnow,])-m)
    Sb     = Sb + Nk*outer(mdiff,mdiff)
  }
  # 3. cost function
  costS = Sb-C*Sw
  # 4. use top
  projection = aux.adjqr(RSpectra::eigs(costS, ndim)$vectors)


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
