#' Linear Discriminant Analysis
#'
#' Linear Discriminant Analysis (LDA) originally aims to find a set of features
#' that best separate groups of data. Since we need \emph{label} information,
#' LDA belongs to a class of supervised methods of performing classification.
#' However, since it is based on finding \emph{suitable} projections, it can still
#' be used to do dimension reduction. We support both binary and multiple-class cases.
#' Note that the target dimension \code{ndim} should be \emph{less than or equal to} \code{K-1},
#' where \code{K} is the number of classes, or \code{K=length(unique(label))}. Our code
#' automatically gives bounds on user's choice to correspond to what theory has shown. See
#' the comments section for more details.
#'
#' @section Limit of Target Dimension Selection:
#' In unsupervised algorithms, selection of \code{ndim} is arbitrary as long as
#' the target dimension is lower-dimensional than original data dimension, i.e., \code{ndim < p}.
#' In LDA, it is \emph{not allowed}. Suppose we have \code{K} classes, then its formulation on
#' \eqn{S_B}, between-group variance, has maximum rank of \code{K-1}. Therefore, the maximal
#' subspace can only be spanned by at most \code{K-1} orthogonal vectors.
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=33)-100
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+100
#'
#' ## merge the data and create a label correspondingly
#' Y      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## perform onto 2-dimensional space
#' output = do.lda(Y, label, ndim=2)
#'
#' ## visualize
#' plot(output$Y[,1], output$Y[,2], main="3 groups on 2d plane")
#'
#' @references
#' \insertRef{fisher_use_1936}{Rdimtools}
#'
#' \insertRef{fukunaga_introduction_1990}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LDA
#' @export
do.lda <- function(X, label, ndim=2){
  ## Preprocessing
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label vector
  if ((!is.vector(label))||(length(label)!=n)){
    stop("* do.lda : 'label' is required to be a vector of class labels.")
  }
  label  = as.numeric(as.factor(label))
  ulabel = unique(label)
  K      = length(ulabel)
  if (K==1){
    stop("* do.lda : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    warning("* do.lda : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }

  #   3. ndim
  if (!check_ndim(ndim,p)){
    stop("* do.lda : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  if (ndim>=K){
    warning("* do.lda : by the nature of LDA, target dimension 'ndim' is adjusted to match maximally permissible subspace.")
    ndim = (K-1)
  }
  #   4. perform CENTERING
  tmplist = aux.preprocess(X,type="center")
  trfinfo = tmplist$info
  trfinfo$algtype = "linear"
  pX      = tmplist$pX

  ## Main Computation
  result = list()
  if (K==2){ ## 2-class
    idx1  = which(label==ulabel[1])
    idx2  = which(label==ulabel[2])
    SW    = lda_outer(pX[idx1,]) + lda_outer(pX[idx2,])
    mdiff = matrix(colMeans(pX[idx2,])-colMeans(pX[idx1,]))
    RLIN  = Rlinsolve::lsolve.bicgstab(SW, mdiff, verbose=FALSE)
    w     = aux.adjprojection(as.matrix(RLIN$x))
    Y     = pX%*%w

    result$Y = Y
    result$trfinfo = trfinfo
    result$projection = w
  } else { ## K-class : maximally, (K-1) dimension possible, you know I'm saying?
    # 1. compute S_W  : within-group variance for multiclss problem
    SW = array(0,c(p,p))
    for (i in 1:K){
      idxnow = which(label==ulabel[i])
      SW     = SW + lda_outer(pX[idxnow,])
    }
    # 2. compute S_B  : between-group variance for multiclass problem
    SB = array(0,c(p,p))
    m  = colMeans(pX)
    for (i in 1:K){
      idxnow = which(label==ulabel[i])
      Nk     = length(idxnow)
      mdiff  = (colMeans(pX[idxnow,])-m)
      SB     = SB + Nk*outer(mdiff,mdiff)
    }
    RLIN = Rlinsolve::lsolve.bicgstab(SW, SB, verbose=FALSE)
    W    = RLIN$x
    topW = aux.adjprojection(RSpectra::eigs(W, ndim)$vectors)
    Y    = pX%*%topW

    result$Y = Y
    result$trfinfo = trfinfo
    result$projection = topW
  }
  ## Return results
  return(result)
}


## sum of outer products to calculate S_B and S_W
#' @keywords internal
#' @noRd
lda_outer <- function(X){
  p = ncol(X)
  output = array(0,c(p,p))
  for (i in 1:nrow(X)){
    output = output + outer(X[i,],X[i,])
  }
  return(output)
}
