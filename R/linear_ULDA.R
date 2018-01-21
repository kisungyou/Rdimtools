#' Uncorrelated Linear Discriminant Analysis
#'
#'
#' @examples
#' \dontrun{
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=33)-100
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+100
#'
#' ## merge the data and create a label correspondingly
#' X      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## compare with LDA
#' out1 = do.lda(X, label)
#' out2 = do.ulda(X, label)
#'
#' ## visualize
#' par(mfrow=c(1,2))
#' plot(out1$Y[,1], out1$Y[,2])
#' plot(out2$Y[,1], out2$Y[,2])
#' }
#'
#' @rdname linear_ULDA
#' @export
do.ulda <- function(X, label, ndim=2, preprocess=c("center","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label vector
  if ((!is.vector(label))||(length(label)!=n)){
    stop("* do.ulda : 'label' is required to be a vector of class labels.")
  }
  label  = as.numeric(as.factor(label))
  ulabel = unique(label)
  K      = length(ulabel)
  if (K==1){
    stop("* do.ulda : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    stop("* do.ulda : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.ulda : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"

  #   2. Sb and Sw
  #   2-1. Sb
  mu_overall = colMeans(pX)
  Sb = array(0,c(p,p))
  for (i in 1:K){
    idxnow = which(label==ulabel[i])
    Pi     = length(idxnow)/n
    mdiff  = (colMeans(pX[idxnow,])-mu_overall)
    Sb     = Sb + Pi*outer(mdiff,mdiff)
  }
  #   2-2. Sw
  Sw = array(0,c(p,p))
  for (i in 1:K){
    idxnow = which(label==ulabel[i])
    Pi     = length(idxnow)/n
    Si     = array(0,c(p,p))
    mu_K   = colMeans(pX[idxnow,])
    for (j in 1:length(idxnow)){
      mdiff = (as.vector(pX[idxnow[j],])-mu_K)
      Si = Si + outer(mdiff,mdiff)
    }
    Sw = Sw + Pi*Si
  }
  #   2-3. pseudo-inverse for Sw
  invSw = pracma::pinv(Sw)

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR UNCORRELATED LDA
  projection = array(0,c(p,ndim))
  #   1. take the first vector
  projection[,1] = as.vector(aux.geigen(Sb, Sw, 1, maximal=TRUE))
  #   2. iterate from 2:ndim
  diagIp = diag(p)
  for (i in 2:ndim){
    # 2-1. select D
    D = t(matrix(projection[,1:(i-1)], nrow=p))
    # 2-2. solve for inverse problem
    tmpA = D%*%invSw%*%t(D)
    tmpB = D%*%invSw
    if (i==2){
      tmpsolve = tmpB/as.double(tmpA)
      Pmat = diagIp - outer(as.vector(D), as.vector(tmpsolve))
    } else {
      tmpsolve = Rlinsolve::lsolve.bicgstab(tmpA, tmpB, verbose=FALSE)$x
      Pmat = diagIp - t(D)%*%tmpsolve
    }
    # 2-4. cost function for outer generalized eigenvalue problem and solve
    projection[,i] = as.vector(aux.geigen(Pmat%*%Sb, Sw, 1, maximal=TRUE))
  }

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
