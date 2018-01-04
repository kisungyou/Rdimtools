#' Modified Orthogonal Discriminant Projection
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center" and two other options "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.
#' @param alpha balancing parameter of non-local and local scatter in \eqn{[0,1]}.
#' @param beta scaling control parameter for distant pairs of data in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' ## generate 3 different groups of data X and label vector
#' x1 = matrix(rnorm(4*10), nrow=10)-20
#' x2 = matrix(rnorm(4*10), nrow=10)
#' x3 = matrix(rnorm(4*10), nrow=10)+20
#' X  = rbind(x1, x2, x3)
#' label = c(rep(1,10), rep(2,10), rep(3,10))
#'
#' ## try different beta (scaling control) parameter
#' out1 = do.modp(X, label, beta=1)
#' out2 = do.modp(X, label, beta=10)
#' out3 = do.modp(X, label, beta=100)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="modp::beta=1")
#' plot(out2$Y[,1], out2$Y[,2], main="modp::beta=10")
#' plot(out3$Y[,1], out3$Y[,2], main="modp::beta=100")
#'
#' @references
#' \insertRef{zhang_modified_2011}{Rdimtools}
#'
#' @rdname linear_MODP
#' @export
do.modp <- function(X, label, ndim=2, preprocess=c("center","decorrelate","whiten"),
                   type=c("proportion",0.1), symmetric=c("union","intersect","asymmetric"),
                   alpha = 0.5, beta = 10){
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
      stop("* do.modp : no degerate class of size 1 is allowed.")
    }
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    warning("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.modp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. nbd-type
  nbdtype = type
  #   6. nbd-symmetric
  if (missing(symmetric)){
    nbdsymmetric = "union"
  } else {
    nbdsymmetric = match.arg(symmetric)
  }
  #   7. alpha and beta
  alpha = as.double(alpha)
  if (!check_NumMM(alpha,0,1,compact=TRUE)){stop("* do.modp : 'alpha' is a balancing parameter in [0,1].")}
  beta = as.double(beta)
  if (!check_NumMM(beta,0,Inf,compact=FALSE)){stop("* do.modp : 'beta' is a scaling control parameter in (0,inf).")}
  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. Preprocessing the data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
  #   2. neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask
  #   3. Distance Matrix Squared
  Dmat2 = (as.matrix(dist(pX))^2)
  #   4. Construct W : weight matrix
  W = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if ((nbdmask[i,j]==TRUE)&&(nbdmask[j,i]==TRUE)){ # neighbors of each other
        if (label[i]==label[j]){
          expval = exp((-Dmat2[i,j])/beta)
          W[i,j] = expval
          W[j,i] = expval
        } else {
          expval = exp((-Dmat2[i,j])/beta) ### this part should be changed with Modified ODP
          comval = expval*abs(cor(as.vector(pX[i,]), as.vector(pX[j,])))
          W[i,j] = comval
          W[j,i] = comval
        }
      }
    }
  }
  #   5. Construct Sl and St
  #   5-1. Sl : local
  L  = diag(rowSums(W))-W
  Sl = (t(pX)%*%L%*%pX)/(2*n*n)
  #   5-2. St : total : non-local Sn = St-Sl
  St = aux_scatter_pairwise(pX)/(2*n*n)
  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN ODP
  #   1. cost function
  costS = ((1-alpha)*St)-(alpha*Sl)
  #   2. top eigenvectors
  projection = RSpectra::eigs(costS, ndim)$vectors
  #   3. adjust projection matrix
  projection = aux.adjprojection(projection)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}





