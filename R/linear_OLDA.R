#' Orthogonal Linear Discriminant Analysis
#'
#' Orthogonal LDA (OLDA) is an extension of classical LDA where the discriminant vectors are
#' orthogonal to each other.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
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
#' X     = as.matrix(iris[,1:4])
#' label = as.integer(iris$Species)
#'
#' ## compare with LDA
#' out1 = do.lda(X, label)
#' out2 = do.olda(X, label)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(out1$Y, col=label, main="LDA")
#' plot(out2$Y, col=label, main="Orthogonal LDA")
#' par(opar)
#'
#' @references
#' \insertRef{ye_characterization_2005}{Rdimtools}
#'
#' @rdname linear_OLDA
#' @author Kisung You
#' @concept linear_methods 
#' @export
do.olda <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label vector
  label  = check_label(label, n)
  ulabel = unique(label)
  K      = length(ulabel)
  if (K==1){
    stop("* do.olda : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    stop("* do.olda : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.olda : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. Sb and Sw
  #   2-1. Sb
  mu_overall = colMeans(pX)
  Sb = array(0,c(p,p))
  for (i in 1:K){
    idxnow = which(label==ulabel[i])
    Pi     = length(idxnow)/n
    mdiff  = (colMeans(pX[idxnow,])-mu_overall)
    Sb     = Sb + Pi*outer(mdiff,mdiff)
  }
  #   2-2. Sw
  Sw = array(0,c(p,p))
  for (i in 1:K){
    idxnow = which(label==ulabel[i])
    Pi     = length(idxnow)/n
    Si     = array(0,c(p,p))
    mu_K   = colMeans(pX[idxnow,])
    for (j in 1:length(idxnow)){
      mdiff = (as.vector(pX[idxnow[j],])-mu_K)
      Si = Si + outer(mdiff,mdiff)
    }
    Sw = Sw + Pi*Si
  }
  #   2-3. pseudo-inverse for Sw; using my function
  invSw = aux.pinv(Sw)

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR UNCORRELATED LDA
  #   1. initialize
  Dmat = as.vector(aux.geigen(Sb, Sw, 1, maximal=TRUE))
  #   2. step-by-step computation
  diagIp = diag(p)
  for (i in 2:ndim){
    #   2-1.
    if (i==2){
      Dmat = matrix(Dmat)
    }
    tmpA = t(Dmat)%*%invSw%*%Dmat
    tmpB = t(Dmat)%*%invSw

    #   2-2. solve intermediate inverse problem
    tmpsolve = aux.bicgstab(tmpA, tmpB, verbose=FALSE)$x
    Pmat = diagIp - (Dmat)%*%tmpsolve

    # 2-3. cost function for outer generalized eigenvalue problem and solve
    csolution = as.vector(aux.geigen(Pmat%*%Sb, Sw, 1, maximal=TRUE))
    Dmat = cbind(Dmat, csolution)
  }
  # 2-4. remove column names
  colnames(Dmat)=NULL

  #------------------------------------------------------------------------
  ## RETURN
  #   1. adjust with orthogonalization
  projection = aux.adjqr(Dmat)

  #   2. return output
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
