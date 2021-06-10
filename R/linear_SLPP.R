#' Supervised Locality Preserving Projection
#'
#' As its names suggests, Supervised Locality Preserving Projection (SLPP) is a variant of LPP
#' in that it replaces neighborhood network construction schematic with class information in that
#' if two nodes belong to the same class, it assigns weight of 1, i.e., \eqn{S_{ij}=1} if \eqn{x_i} and
#' \eqn{x_j} have same class labelings.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center" and other options of "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150, 50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' ## compare SLPP with LPP
#' outLPP  <- do.lpp(X)
#' outSLPP <- do.slpp(X, label)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(outLPP$Y,  pch=19, col=label, main="LPP")
#' plot(outSLPP$Y, pch=19, col=label, main="SLPP")
#' par(opar)
#'
#' @references
#' \insertRef{zheng_gabor_2007}{Rdimtools}
#'
#' @seealso \code{\link{do.lpp}}
#' @author Kisung You
#' @rdname linear_SLPP
#' @concept linear_methods
#' @export
do.slpp <- function(X, label, ndim=2, preprocess=c("center","decorrelate","whiten")){
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
      stop("* do.slpp : no degerate class of size 1 is allowed.")
    }
  }
  if (any(is.na(label))||(any(is.infinite(label)))){stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){    stop("* do.slpp : 'ndim' is a positive integer in [1,#(covariates)].")  }
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  } else {    algpreprocess = match.arg(preprocess)  }
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"

  #------------------------------------------------------------------------
  ## COMPUTATION PART
  #   1. PCA Preprocessing
  eigtest = eigen(cov(pX), only.values=TRUE)
  pcadim  = sum(eigtest$values > 0)

  if (pcadim <= ndim){
    warning("* do.slpp : target 'ndim' is larger than intrinsic data dimension achieved from PCA.")
    projection_first = diag(p)
    pcapX = pX
  } else{
    projection_first = aux.adjprojection(eigen(cov(pX))$vectors[,1:pcadim])
    pcapX = pX%*%projection_first
  }

  #   2. computation following without ordering
  S = array(0,c(n,n))
  for (i in 1:length(ulabel)){
    # 2-1. find the target label and matching index
    tgtidx = which(label==ulabel[i])
    # 2-2. generate allones-{but diagonal}
    tmpslppmat = slpp_allones(length(tgtidx))
    # 2-3. plug-in
    S[tgtidx,tgtidx] = tmpslppmat
  }
  diag(S)=0

  #   3. main LPP unit
  D = diag(rowSums(S))
  L = D-S

  LHS = t(pcapX)%*%L%*%pcapX
  RHS = t(pcapX)%*%D%*%pcapX

  #   4. use lowest vectors
  projection_second = aux.geigen(LHS, RHS, ndim, maximal=FALSE)


  #------------------------------------------------------------------------
  ## RETURN
  #   1. adjust projection
  projection = (projection_first%*%projection_second)
  projection = aux.adjprojection(projection)

  #   2.
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}




# : -----------------------------------------------------------------------
#' @keywords internal
#' @noRd
slpp_allones <- function(n){
  output = matrix(rep(1,n^2), nrow=n)
  diag(output) = 0
  return(output)
}
