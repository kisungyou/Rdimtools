#' Local Fisher Discriminant Analysis
#'
#' Local Fisher Discriminant Analysis (LFDA) is a linear dimension reduction method for
#' supervised case, i.e., labels are given. It reflects \emph{local} information to overcome
#' undesired results of traditional Fisher Discriminant Analysis which results in a poor mapping
#' when samples in a single class form form several separate clusters.
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.
#' @param localscaling \code{TRUE} to use local scaling method for construction affinity matrix, \code{FALSE} for binary affinity.
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
#' X     = rbind(x1, x2, x3)
#' label = rep(1:3, each=10)
#'
#' ## try different affinity matrices
#' out1 = do.lfda(X, label)
#' out2 = do.lfda(X, label, localscaling=FALSE)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(out1$Y, col=label, main="binary affinity matrix")
#' plot(out2$Y, col=label, main="local scaling affinity")
#' par(opar)
#'
#' @references
#' \insertRef{sugiyama_local_2006}{Rdimtools}
#'
#' \insertRef{zelnik-manor_selftuning_2005}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LFDA
#' @concept linear_methods
#' @export
do.lfda <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"),
                    type=c("proportion",0.1), symmetric=c("union","intersect","asymmetric"),
                    localscaling=TRUE){
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
      stop("* do.lfda : no degerate class of size 1 is allowed.")
    }
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lfda : 'ndim' is a positive integer in [1,#(covariates)).")}
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
  #   7. localscaling
  if (!is.logical(localscaling)){
    stop("* do.lfda : 'localscaling' must be a logical flag.")
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. Preprocessing the data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #   2. neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask
  #   3. construct A : neighborhood matrix
  #      localscaling / binary affinity
  if (localscaling==FALSE){
    A = nbdmask * 1.0
  } else {
    #   3-1. compute sigma_i
    vec_sigmai = rep(0,n)
    for (i in 1:n){
      #   3-2. select same class
      tgtidxAi = which(nbdmask[i,])
      #   3-3. compute
      if (length(tgtidxAi)<1){
        message(paste("* do.lfda : ",i,"-th element has no neighbors."))
        vec_sigmai[i] = 0.001;
      } else if (length(tgtidxAi)==1){
        vecdiff1 = as.vector(pX[i,])-as.vector(pX[tgtidxAi,])
        vec_sigmai[i] = sqrt(sum(vecdiff1*vecdiff1))
      } else {
        simplevec = as.vector(pX[i,])
        simplemat = as.matrix(pX[tgtidxAi,])
        vec_sigmai[i] = method_lfda_maximaldistance(simplevec, simplemat)
      }
    }
    Dmat = (as.matrix(dist(pX))^2)
    A = exp(-diag(1/vec_sigmai)%*%Dmat%*%diag(1/vec_sigmai))
  }
  #   4. class-wise elements counting
  ulabel_counts = rep(0,length(ulabel))
  ulabel_idx    = list()
  for (i in 1:length(ulabel)){
    whichlabel = which(label==ulabel[i])
    ulabel_counts[i] = length(whichlabel)
    ulabel_idx[[i]]  = whichlabel
  }
  #   5. Construct Aijw and Aijb
  Aijw = array(0,c(n,n))
  Aijb = array(0,c(n,n))
  for (i in 1:n){
    ylab1 = label[i]
    for (j in 1:n){
      ylab2 = label[j]
      if (ylab1==ylab2){
        nc = ulabel_counts[which(ulabel==ylab1)]
        Aijw[i,j] = A[i,j]/nc
        Aijb[i,j] = A[i,j]*((1/n)-(1/nc))
      } else {
        Aijb[i,j] = 1/n
      }
    }
  }
  #   6. Construct Swbar and Sbbar
  Swbar = array(0,c(p,p))
  Sbbar = array(0,c(p,p))
  for (i in 1:n){
    vec1 = as.vector(pX[i,])
    for (j in 1:n){
      vec2 = as.vector(pX[j,])
      vecdiff = vec1-vec2
      outdiff = outer(vecdiff,vecdiff)

      Swbar = Swbar + ((Aijw[i,j])*outdiff)
      Sbbar = Sbbar + ((Aijb[i,j])*outdiff)
    }
  }


  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN LFDA : USE TOP EIGENVECTORS
  Smbar = Swbar+Sbbar
  projection = aux.geigen(Smbar, Swbar, ndim, maximal=TRUE)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
