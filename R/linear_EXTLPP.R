#' Extended Locality Preserving Projection
#'
#'
#'
#' @references
#' \insertRef{shikkenawis_improving_2012}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_EXTLPP
#' @export
do.extlpp <- function(X, ndim=2, numk=max(ceiling(nrow(X)/10),2),
                      preprocess=c("center","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.extlpp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. numk
  numk = as.integer(numk)
  if (!check_NumMM(numk,1,n/2,compact=FALSE)){stop("* do.extlpp : 'numk' should be an integer in [2,nrow(X)/2).")}
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. preprocessing
  if (algpreprocess=="null"){
    trfinfo = list()
    trfinfo$type = "null"
    pX = as.matrix(X,nrow=nrow(X));
  } else {
    tmplist = aux.preprocess(X,type=algpreprocess)
    trfinfo = tmplist$info
    pX      = tmplist$pX
  }
  trfinfo$algtype = "linear"

  #   2. K-Means Clustering
  kclust     = stats::kmeans(pX, numk)
  clustlabel = kclust$cluster
  clustidx   = list() # for each label, find the corresponding # length-'numk' list
  for (i in 1:numk){
    clustidx[[i]] = which(clustlabel==unique(clustlabel)[i])
  }

  #   3. pairwise distance
  PD = as.matrix(dist(pX))
  vecb = rep(0,numk)
  for (i in 1:numk){
    tgtidx = clustidx[[i]]
    vecb[i] = max(PD[tgtidx,tgtidx])
  }
  veca = rep(min(vecb)/20,numk)

  #   4. compute S
  S = array(0,c(n,n))
  for (i in 1:numk){
    tgtidx = clustidx[[i]]
    S[tgtidx,tgtidx] = method_trfextlpp(PD[tgtidx,tgtidx],veca[i],vecb[i])
  }
  diag(S) = 0.0

  #   5. graph laplaciana and generalized eigenvalue problem
  D = diag(rowSums(S))
  L = D-S

  LHS = t(pX)%*%L%*%pX
  RHS = t(pX)%*%D%*%pX

  #   6. compute Projection Matrix
  geigs = geigen::geigen(LHS, RHS, TRUE)
  projection = matrix(geigs$vectors[,1:ndim],nrow=p)
  eigenvalue = as.vector(geigs$values[1:ndim])

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$eigval = eigenvalue
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}

