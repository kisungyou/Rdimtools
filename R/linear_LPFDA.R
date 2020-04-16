#' Locality Preserving Fisher Discriminant Analysis
#'
#' Locality Preserving Fisher Discriminant Analysis (LPFDA) is a supervised variant of LPP.
#' It can also be seemed as an improved version of LDA where the locality structure of the data
#' is preserved. The algorithm aims at getting a subspace projection matrix by solving a generalized
#' eigenvalue problem.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param t bandwidth parameter for heat kernel in \eqn{(0,\infty)}.
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
#' ## try different proportion of connected edges
#' out1 = do.lpfda(X, label, type=c("proportion",0.01))
#' out2 = do.lpfda(X, label, type=c("proportion",0.1))
#' out3 = do.lpfda(X, label, type=c("proportion",0.25))
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="1% connectivity")
#' plot(out2$Y, main="10% connectivity")
#' plot(out3$Y, main="25% connectivity")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zhao_locality_2009}{Rdimtools}
#'
#' @rdname linear_LPFDA
#' @author Kisung You
#' @concept linear_methods 
#' @export
do.lpfda <- function(X, label, ndim=2, type=c("proportion",0.1),
                     preprocess=c("center","scale","cscale","whiten","decorrelate"), t=10.0){
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
    stop("* do.lpfda : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    stop("* do.lpfda : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lpfda : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. type
  nbdtype = type
  nbdsymmetric = "union"
  #   5. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   6. t : kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t, .Machine$double.eps*10, Inf)){stop("* do.lpfda : 't' is a positive kernel bandwidth parameter.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  # 2. process : neighborhood selection
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR LPFDA
  # 1. LPP part : compute weight matrix W and other auxiliary's
  Dsqmat = exp(-(as.matrix(dist(pX))^2)/t)
  W = Dsqmat*nbdmask
  L = diag(rowSums(W))-W
  # 2. FDA part : compute two scatter matrices
  # 2-1. compute S_W  : within-group variance for multiclss problem
  Sw = array(0,c(p,p))
  for (i in 1:K){
    idxnow = which(label==ulabel[i])
    Sw     = Sw + lda_outer(pX[idxnow,])
  }
  # 2-2. compute S_B  : between-group variance for multiclass problem
  Sb = array(0,c(p,p))
  m  = colMeans(pX)
  for (i in 1:K){
    idxnow = which(label==ulabel[i])
    Nk     = length(idxnow)
    mdiff  = (colMeans(pX[idxnow,])-m)
    Sb     = Sb + Nk*outer(mdiff,mdiff)
  }
  # 3. cost function for geigen : use Maximal
  costL      = Sb - t(pX)%*%L%*%pX
  projection = aux.geigen(costL, Sw, ndim, maximal=TRUE)


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
