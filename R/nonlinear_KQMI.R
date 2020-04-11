#' Kernel Quadratic Mutual Information
#'
#' Kernel Quadratic Mutual Information (KQMI) is a supervised linear dimension reduction method.
#' Quadratic Mutual Information is an efficient nonparametric estimation method for Mutual Information
#' for class labels not requiring class priors. The method re-states the estimation procedure in terms of
#' kernel objective in the graph embedding framework.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
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
#' ## generate 3 different groups of data X and label vector
#' x1 = matrix(rnorm(4*10), nrow=10)-20
#' x2 = matrix(rnorm(4*10), nrow=10)
#' x3 = matrix(rnorm(4*10), nrow=10)+20
#' X  = rbind(x1, x2, x3)
#' label = c(rep(1,10), rep(2,10), rep(3,10))
#'
#' ## try different kernel bandwidths
#' out1 = do.kqmi(X, label, t=0.01)
#' out2 = do.kqmi(X, label, t=1)
#' out3 = do.kqmi(X, label, t=100)
#'
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, main="KQMI::t=0.01")
#' plot(out2$Y, col=label, main="KQMI::t=1")
#' plot(out3$Y, col=label, main="KQMI::t=100")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{bouzas_graph_2015}{Rdimtools}
#'
#' @seealso \code{\link{do.lqmi}}
#' @author Kisung You
#' @rdname nonlinear_KQMI
#' @export
do.kqmi <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","whiten","decorrelate"), t=10){
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
    message("* do.kqmi : the method is intended for case when ndim=(C-1).")
  }
  if (nlabel==1){
    stop("* do.kqmi : 'label' should have at least 2 unique labelings.")
  }
  if (nlabel==n){
    stop("* do.kqmi : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.kqmi : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. t : kernel bandwidth parameter
  t = as.double(t)
  if (!check_NumMM(t, .Machine$double.eps, Inf, compact=FALSE)){stop("* do.kqmi : 't' is a bandwidth parameter for gaussian kernel.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN STEPS FOR Kernel QMI
  #   1. compute M and symmetrize
  M    = qmi_M(pX, label)
  symM = ((M+t(M))/2)
  #   2. kernel and centered version
  Ktilde = exp(-(as.matrix(dist(pX))^2)/t)
  if ((any(is.na(Ktilde)))||(any(is.infinite(Ktilde)))){
    stop("* do.kqmi : kernel matrix contains Inf or NA. Please adjust kernel bandwidth value.")
  }
  En     = array(1/n, c(n,n))
  matK   = (Ktilde - (En%*%Ktilde) - (Ktilde%*%En) + (En%*%Ktilde%*%En))
  #   3. decompose kernel matrix
  eigmatK = base::eigen(matK)
  vecL    = eigmatK$values
  matP    = eigmatK$vectors
  #   4. solve the eigenvalue problem
  costmat = (t(matP)%*%symM%*%matP)
  eigcost = base::eigen(costmat)
  matB    = eigcost$vectors
  #   5. compute optimal projection
  invvecL = 1/vecL
  nogoodL = which(is.na(invvecL)||is.infinite(invvecL))
  if (length(nogoodL)>0){
    invvecL[nogoodL] = 0
  }
  invmatL = diag(invvecL)
  matA    = (matP%*%invmatL%*%matB)
  #   6. compute pseudo-projection Astar
  pseudoproj = qr.Q(qr(matA))[,1:ndim]


  #------------------------------------------------------------------------
  ## RETURN THE RESULTS
  result = list()
  result$Y = matK%*%pseudoproj
  result$trfinfo = trfinfo
  return(result)
}
