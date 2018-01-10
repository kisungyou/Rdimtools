#' Localized Sliced Inverse Regression
#'
#'
#'
#' @seealso \code{\link{do.sir}}
#' @author Kisung You
#' @rdname linear_LSIR
#' @export
do.lsir <- function(X, response, ndim=2, h=max(2, round(nrow(X)/5)), preprocess=c("center","decorrelate","whiten"),
                    numk = max(2, round(nrow(X)/10)), tau=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. response
  response = as.double(response)
  if ((!is.vector(response))||(any(is.na(response)))){
    stop("* do.sir : 'response' should be a vector containing no NA values.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lsir : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. h : number of slices
  h = as.integer(h)
  if (!check_NumMM(h,2,ceiling(n/2),compact=TRUE)){stop("* do.lsir : the number of slices should be in [2,n/2].")}
  #   5. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   6. numk
  numk = as.integer(numk)
  if (!check_NumMM(numk,1,(n-1),compact=TRUE)){stop("* do.lsir : 'numk' should be in [1,n-1].")}
  #   7. tau
  tau = as.double(tau)
  if (!check_NumMM(tau,0,Inf,compact=TRUE)){stop("* do.lsir : 'tau' should be a nonnegative real number.")}


  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
  #   2. build label matrix
  label  = as.integer(sir_makelabel(response, h))
  ulabel = unique(label)
  nlabel = length(ulabel)
  #   3. compute classwise and overall mean
  class_mean  = array(0,c(nlabel,p))
  class_count = rep(0,nlabel)
  for (i in 1:nlabel){
    idxclass = which(label==ulabel[i])
    class_mean[i,] = as.vector(colMeans(pX[idxclass,]) )
    class_count[i] = length(idxclass)
  }
  all_mean = as.vector(colMeans(pX))
  #   4. compute Empirical Covariance
  mat_Sigma = aux_scatter(pX, all_mean)/n
  #   5. compute Empirical Between-ness Covariance with localization
  mean_local = array(0,c(n,p))
  for (i in 1:n){
    # 5-1. target label index
    tgtlabel = setdiff(which(label==label[i]), i)
    # 5-2. get the logical of smallest distances
    tvec = as.vector(pX[i,])
    tmat = pX[tgtlabel,]
    partidx  = lsir_smallest(tvec, tmat, numk)
    partlabel = tgtlabel[partidx]
    # 5-3. compute the mean
    mean_local[i,] = as.vector(colMeans(pX[partlabel,]))
  }
  # 5-4. compute mean scatter
  mat_Gamma = aux_scatter(mean_local, all_mean)/n + tau*diag(p)


  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION
  #   1. do matrix inversion.. I hate it.
  costInv = lsolve.bicgstab(mat_Sigma, mat_Gamma, verbose=FALSE)$x
  #   2. find top eigenvectors
  projection = aux.adjprojection(RSpectra::eigs(costInv, ndim)$vectors)


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}




#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
lsir_smallest <- function(tvec, tmat, numk){
  n = nrow(tmat)
  distvec = rep(0,n)
  for (i in 1:n){
    vecdiff = tvec-as.vector(tmat[i,])
    distvec[i] = sum(vecdiff*vecdiff)
  }
  orderdist = order(distvec)
  orderpart = orderdist[1:max(min(n,numk),1)]

  logicalvec = rep(FALSE,n)
  logicalvec[orderpart] = TRUE
  return(logicalvec)
}
