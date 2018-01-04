#' Kernel Local Fisher Discriminant Analysis
#'
#' Kernel LFDA is a nonlinear extension of LFDA method using kernel trick. It applies conventional kernel method
#' to extend excavation of hidden patterns in a more flexible manner in tradeoff of computational load. For simplicity,
#' only the gaussian kernel parametrized by its bandwidth \code{t} is supported.
#'
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
#' @param localscaling \code{TRUE} to use local scaling method for construction affinity matrix, \code{FALSE} for binary affinity.
#' @param t bandwidth parameter for heat kernel in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate 3 different groups of data X and label vector
#' x1 = matrix(rnorm(4*10), nrow=10)-20
#' x2 = matrix(rnorm(4*10), nrow=10)
#' x3 = matrix(rnorm(4*10), nrow=10)+20
#' X  = rbind(x1, x2, x3)
#' label = c(rep(1,10), rep(2,10), rep(3,10))
#'
#' ## try different affinity matrices
#' out1 = do.klfda(X, label, t=0.1)
#' out2 = do.klfda(X, label, t=1)
#' out3 = do.klfda(X, label, t=10)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="bandwidth=0.1")
#' plot(out2$Y[,1], out2$Y[,2], main="bandwidth=1")
#' plot(out3$Y[,1], out3$Y[,2], main="bandwidth=10")
#' }
#'
#' @references
#' \insertRef{sugiyama_local_2006}{Rdimtools}
#'
#' \insertRef{zelnik-manor_self-tuning_2005}{Rdimtools}
#'
#' @seealso \code{\link{do.lfda}}
#' @author Kisung You
#' @rdname nonlinear_KLFDA
#' @export
do.klfda <- function(X, label, ndim=2, preprocess=c("center","decorrelate","whiten"),
                     type=c("proportion",0.1), symmetric=c("union","intersect","asymmetric"),
                     localscaling=TRUE, t=1.0){
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
      stop("* do.klfda : no degerate class of size 1 is allowed.")
    }
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    warning("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.klfda : 'ndim' is a positive integer in [1,#(covariates)).")}
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
  #   8. t : kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t, 0, 1e+10, compact=FALSE)){stop("* do.klfda : 't' is a bandwidth parameter for gaussian kernel.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. Preprocessing the data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "nonlinear"
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
  Aijm = Aijw + Aijb
  #   6. Compute Laplacian matrices Lijm and Lijb
  Lijw = diag(rowSums(Aijw))-Aijw
  Lijm = diag(rowSums(Aijm))-Aijm
  #   7. compute Kernel Matrix
  K = exp(-(as.matrix(dist(pX))^2)/(2*(t^2)))

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN KLFDA
  #  since direct computation is difficult, I used a detour using Rlinsolve and RSpectra
  LHS = K%*%Lijm%*%K
  RHS = K%*%Lijw%*%K

  CHS = Rlinsolve::lsolve.bicgstab(RHS, LHS, verbose=FALSE)$x
  pseudoproj = RSpectra::eigs(CHS, ndim)$vectors

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = K%*%pseudoproj
  result$trfinfo = trfinfo
  return(result)
}
