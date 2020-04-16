#' Regularized Linear Discriminant Analysis
#'
#' In small sample case, Linear Discriminant Analysis (LDA) may suffer from
#' rank deficiency issue. Applied mathematics has used Tikhonov regularization -
#' also known as \eqn{\ell_2} regularization/shrinkage - to adjust linear operator.
#' Regularized Linear Discriminant Analysis (RLDA) adopts such idea to stabilize
#' eigendecomposition in LDA formulation.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param alpha Tikhonow regularization parameter.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \dontrun{
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.factor(iris$Species)
#'
#' ## try different regularization parameters
#' out1 <- do.rlda(X, label, alpha=0.001)
#' out2 <- do.rlda(X, label, alpha=0.01)
#' out3 <- do.rlda(X, label, alpha=100)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, main="RLDA::alpha=0.1")
#' plot(out2$Y, col=label, main="RLDA::alpha=1")
#' plot(out3$Y, col=label, main="RLDA::alpha=10")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{friedman_regularized_1989}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_RLDA
#' @export
do.rlda <- function(X, label, ndim=2, alpha=1.0){
  ## Note : refer to do.klfda
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
      stop("* do.rlda : no degerate class of size 1 is allowed.")
    }
  }
  K      = length(ulabel)
  if (K==1){
    stop("* do.rlda : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    warning("* do.rlda : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.rlda : 'ndim' is a positive integer in [1,#(covariates)).")}

  #   4. alpha : regularization parameter
  alpha = as.double(alpha)
  if (alpha==0){
    stop("* do.rlda : 'alpha=0' condition leads to applying lda. Use 'do.lda' instead.")
  }
  if (!check_NumMM(alpha,0,Inf,compact=TRUE)){stop("* do.rlda : 'alpha' is a regularization parameter in (0,Inf).")}

  #   (implicit) : preprocess
  algpreprocess = "center"

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. Preprocessing the data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
  #   2. Sb and St : instead of Sw, we can use St on the denominator.
  #   2-1. St
  mu_overall = as.vector(colMeans(pX))
  St = aux_scatter(pX, mu_overall)
  #   2-2. Sb
  Sb = array(0,c(p,p))
  m  = colMeans(pX)
  for (i in 1:K){
    idxnow = which(label==ulabel[i])
    Nk     = length(idxnow)
    mdiff  = (colMeans(pX[idxnow,])-m)
    Sb     = Sb + Nk*outer(mdiff,mdiff)
  }
  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN rlda
  #  let's use Rlinsolve and geigen structure
  LHS = Sb
  RHS = St + alpha*diag(p)
  #   run Rlinsolve
  W = aux.bicgstab(RHS, LHS, verbose=FALSE)$x
  #   adjust
  topW = aux.adjprojection(RSpectra::eigs(W, ndim)$vectors)
  topW = matrix(as.double(topW), nrow=p)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%topW
  result$trfinfo = trfinfo
  result$projection = topW
  return(result)
}
