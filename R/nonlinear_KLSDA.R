#' Kernel Locality Sensitive Discriminant Analysis
#'
#' Kernel LSDA (KLSDA) is a nonlinear extension of LFDA method using kernel trick. It applies conventional kernel method
#' to extend excavation of hidden patterns in a more flexible manner in tradeoff of computational load. For simplicity,
#' only the gaussian kernel parametrized by its bandwidth \code{t} is supported.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param alpha balancing parameter for between- and within-class scatter in \eqn{[0,1]}.
#' @param k1 the number of same-class neighboring points (homogeneous neighbors).
#' @param k2 the number of different-class neighboring points (heterogeneous neighbors).
#' @param t bandwidth parameter for heat kernel in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' ## generate 3 different groups of data X and label vector
#' x1 = matrix(rnorm(4*10), nrow=10)-50
#' x2 = matrix(rnorm(4*10), nrow=10)
#' x3 = matrix(rnorm(4*10), nrow=10)+50
#' X  = rbind(x1, x2, x3)
#' label = c(rep(1,10), rep(2,10), rep(3,10))
#'
#' ## try different kernel bandwidths
#' out1 = do.klsda(X, label, k1=10, k2=10, t=1)
#' out2 = do.klsda(X, label, k1=10, k2=10, t=5)
#' out3 = do.klsda(X, label, k1=10, k2=10, t=10)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, main="bandwidth=1")
#' plot(out2$Y, col=label, main="bandwidth=15")
#' plot(out3$Y, col=label, main="bandwidth=10")
#' par(opar)
#'
#' @references
#' \insertRef{cai_locality_2007}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_KLSDA
#' @concept nonlinear_methods 
#' @export
do.klsda <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","whiten","decorrelate"),
                    alpha=0.5, k1=max(ceiling(nrow(X)/10),2), k2=max(ceiling(nrow(X)/10),2), t=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label information
  label  = as.numeric(as.factor(label))
  ulabel = unique(label)
  K      = length(ulabel)
  if (K==1){
    stop("* do.klsda : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    stop("* do.klsda : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.klsda : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. alpha
  alpha = as.double(alpha)
  if (!check_NumMM(alpha,0,1)){stop("* do.klsda : 'alpha' is a balancing parameter in [0,1].")}
  #   6. k1 and k2
  k1 = as.integer(k1)
  k2 = as.integer(k2)
  if (!check_NumMM(k1,1,n/2,compact=FALSE)){stop("* do.klsda : 'k1' should be an integer in [2,nrow(X)/2).")}
  if (!check_NumMM(k2,1,n/2,compact=FALSE)){stop("* do.klsda : 'k2' should be an integer in [2,nrow(X)/2).")}
  #   7. t : kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t, 0, 1e+10, compact=FALSE)){stop("* do.klsda : 't' is a bandwidth parameter for gaussian kernel.")}


  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. compute homogeneous (intraclass) and heterogeneous (interclass) neighborhood
  logicalmat = aux.nbdlogical(pX, label, k1, k2)
  Nw = logicalmat$hom
  Nb = logicalmat$het

  #   3. compute Kernel Matrix
  K = exp(-(as.matrix(dist(pX))^2)/(2*(t^2)))

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION FOR LSDA
  #   1. compute auxiliary matrices
  Ww = array(as.logical(Nw+t(Nw)),c(n,n))*1.0; diag(Ww)=0;
  Wb = array(as.logical(Nb+t(Nb)),c(n,n))*1.0; diag(Wb)=0;
  Dw = diag(rowSums(Ww))
  Lb = (diag(rowSums(Wb)) - Wb)

  #   2. make cost function
  LHS = K%*%(alpha*Lb + (1-alpha)*Ww)%*%K
  RHS = K%*%Dw%*%K

  #   3. pseudo-projection; use top eigenvectors
  pseudoproj = aux.geigen(LHS, RHS, ndim, maximal=TRUE)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = K%*%pseudoproj
  result$trfinfo = trfinfo
  return(result)
}
