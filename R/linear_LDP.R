#' Locally Discriminating Projection
#'
#' Locally Discriminating Projection (LDP) is a supervised linear dimension reduction method.
#' It utilizes both label/class information and local neighborhood information to discover
#' the intrinsic structure of the data. It can be considered as an extension
#' of LPP in a supervised manner.
#'
#' @examples
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=33)-100
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+100
#'
#' ## merge the data and create a label correspondingly
#' X      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## try different neighborhood sizes
#' out1 = do.ldp(X, label, type=c("proportion",0.01))
#' out2 = do.ldp(X, label, type=c("proportion",0.05))
#' out3 = do.ldp(X, label, type=c("proportion",0.10))
#'
#' ## visualize
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(out1$Y, col=label, main="1% connectivity")
#' plot(out2$Y, col=label, main="5% connectivity")
#' plot(out3$Y, col=label, main="10% connectivity")
#' par(opar)
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
#' @param beta bandwidth parameter for heat kernel in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @references
#' \insertRef{zhao_local_2006}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LDP
#' @export
do.ldp <- function(X, label, ndim=2, type=c("proportion",0.1),
                   preprocess=c("center","scale","cscale","decorrelate","whiten"), beta=10.0){
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
    stop("* do.ldp : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    stop("* do.ldp : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.ldp : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. type
  nbdtype = type
  nbdsymmetric = "union"
  #   5. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   6. beta
  beta = as.double(beta)
  if (!check_NumMM(beta,0,Inf,compact=FALSE)){stop("* do.ldp : 'beta' is a bandwidth parameter in (0,Inf).")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. build neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN FOR LOCALLY DISCRIMINATING PROJECTION
  #   1. build W
  Dsqmat = exp(-(as.matrix(dist(pX))^2)/beta)
  W = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      alpha = Dsqmat[i,j]
      if (nbdmask[i,j]==TRUE){
        if (label[i]==label[j]){
          thevalue = alpha*(1+alpha)
          W[i,j] = thevalue
          W[j,i] = thevalue
        } else {
          thevalue = alpha*(1-alpha)
          W[i,j] = thevalue
          W[j,i] = thevalue
        }
      }
    }
  }
  #   2. D and L
  D = diag(rowSums(W))
  L = D-W
  #   3. cost function and geigen, BOTTOM solutions
  LHS = t(pX)%*%L%*%pX
  RHS = t(pX)%*%D%*%pX
  projection = aux.geigen(LHS, RHS, ndim, maximal=FALSE)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
