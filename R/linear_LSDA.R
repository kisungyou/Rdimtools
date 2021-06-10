#' Locality Sensitive Discriminant Analysis
#'
#' Locality Sensitive Discriminant Analysis (LSDA) is a supervised linear method.
#' It aims at finding a projection which maximizes the margin between data points from different classes
#' at each local area in which the nearby points with the same label are close to each other while
#' the nearby points with different labels are far apart.
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
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## create a data matrix with clear difference
#' x1 = matrix(rnorm(4*10), nrow=10)-20
#' x2 = matrix(rnorm(4*10), nrow=10)
#' x3 = matrix(rnorm(4*10), nrow=10)+20
#' X  = rbind(x1, x2, x3)
#' label = c(rep(1,10), rep(2,10), rep(3,10))
#'
#' ## try different affinity matrices
#' out1 = do.lsda(X, label, k1=2, k2=2)
#' out2 = do.lsda(X, label, k1=5, k2=5)
#' out3 = do.lsda(X, label, k1=10, k2=10)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, main="nbd size 2")
#' plot(out2$Y, col=label, main="nbd size 5")
#' plot(out3$Y, col=label, main="nbd size 10")
#' par(opar)
#'
#' @references
#' \insertRef{cai_locality_2007}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LSDA
#' @concept linear_methods 
#' @export
do.lsda <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","whiten","decorrelate"),
                    alpha=0.5, k1=max(ceiling(nrow(X)/10),2), k2=max(ceiling(nrow(X)/10),2)){
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
    stop("* do.lsda : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    stop("* do.lsda : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.lsda : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. alpha
  alpha = as.double(alpha)
  if (!check_NumMM(alpha,0,1)){stop("* do.lsda : 'alpha' is a balancing parameter in [0,1].")}
  #   6. k1 and k2
  k1 = as.integer(k1)
  k2 = as.integer(k2)
  if (!check_NumMM(k1,1,n/2,compact=FALSE)){stop("* do.lsda : 'k1' should be an integer in [2,nrow(X)/2).")}
  if (!check_NumMM(k2,1,n/2,compact=FALSE)){stop("* do.lsda : 'k2' should be an integer in [2,nrow(X)/2).")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. compute homogeneous (intraclass) and heterogeneous (interclass) neighborhood
  logicalmat = aux.nbdlogical(pX, label, k1, k2)
  Nw = logicalmat$hom
  Nb = logicalmat$het

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION FOR LSDA
  #   1. compute auxiliary matrices
  Ww = array(as.logical(Nw+t(Nw)),c(n,n))*1.0; diag(Ww)=0;
  Wb = array(as.logical(Nb+t(Nb)),c(n,n))*1.0; diag(Wb)=0;
  Dw = diag(rowSums(Ww))
  Lb = (diag(rowSums(Wb)) - Wb)

  #   2. make cost function
  LHS = t(pX)%*%(alpha*Lb + (1-alpha)*Ww)%*%pX
  RHS = t(pX)%*%Dw%*%pX

  #   3. projection; use top eigenvectors
  projection = aux.geigen(LHS, RHS, ndim, maximal=TRUE)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
