#' Linear Quadratic Mutual Information
#'
#' Linear Quadratic Mutual Information (LQMI) is a supervised linear dimension reduction method.
#' Quadratic Mutual Information is an efficient nonparametric estimation method for Mutual Information
#' for class labels not requiring class priors. For the KQMI formulation, LQMI is a linear equivalent.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
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
#' ## compare against LDA
#' out1 = do.lda(X, label)
#' out2 = do.lqmi(X, label)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(out1$Y, col=label, main="LDA projection")
#' plot(out2$Y, col=label, main="LQMI projection")
#' par(opar)
#'
#' @references
#' \insertRef{bouzas_graph_2015}{Rdimtools}
#'
#' @seealso \code{\link{do.kqmi}}
#' @author Kisung You
#' @rdname linear_LQMI
#' @export
do.lqmi  <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label information
  label  = as.numeric(as.factor(label))
  ulabel = unique(label)
  nlabel = length(ulabel)
  if ((ndim+1)!=nlabel){
    message("* do.lqmi : the method is intended for case when ndim=(C-1).")
  }
  if (nlabel==1){
    stop("* do.lqmi : 'label' should have at least 2 unique labelings.")
  }
  if (nlabel==n){
    stop("* do.lqmi : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.lqmi : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN STEPS FOR LQMI
  #   1. compute M and symmetrize
  M    = qmi_M(pX, label)
  symM = ((M+t(M))/2)
  #   2. eigendecomposition for achieving W
  termA  = t(pX)%*%pX
  termB  = t(pX)%*%symM%*%pX
  linsol = aux.bicgstab(termA, termB, verbose=FALSE)$x
  W      = base::eigen(linsol)$vectors
  #   3. compute projection
  projection = qr.Q(qr(W))[,1:ndim]

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}




# auxiliary for LQMI and KQMI ---------------------------------------------
#' @keywords internal
#' @noRd
qmi_M <- function(X, label){
  n = nrow(X)
  ulabel = unique(label)
  C      = length(ulabel)
  onesN  = rep(1,n)

  # constants
  C_all = 0
  C_btw = rep(0,C)
  for (i in 1:C){
    idxlength = as.double(length(which(label==ulabel[i])))
    C_all     = C_all + (idxlength^2)
    C_btw[i]  = (idxlength/(n^3))
  }
  C_all = C_all/(n^4)
  C_in  = (1/(n^2))

  # matrix
  M = C_all*outer(onesN,onesN)
  for (i in 1:C){
    onesC = (as.vector((label==ulabel[i])*1.0))
    M = M + (C_in*outer(onesC,onesC)) -(2*outer(onesN,as.double(C_btw[i])*onesC))
  }
  return(M)
}
