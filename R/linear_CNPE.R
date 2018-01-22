#' Complete Neighborhood Preserving Embedding
#'
#' @rdname linear_CNPE
#' @author Kisung You
#' @export
do.cnpe <- function(X, ndim=2, type=c("proportion",0.1), preprocess=c("center","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.cnpe : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   3. type
  nbdtype = type
  nbdsymmetric = "union"
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY and LLE step
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"

  #   2. neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)

  #   3. LLE computation
  regparam = 1.0
  W = array(0,c(n,n))
  for (i in 1:n){
    #   3-1. separate target mask vector
    tgtidx  = which(nbdstruct$mask[i,])
    #   3-2. select data
    #        For convenience, target matrix is transposed for Armadillo
    vec_tgt = pX[i,]
    mat_tgt = t(pX[tgtidx,])
    k = ncol(mat_tgt)
    #   3-3. no automatic regularization
    W[i,tgtidx] = method_lleW(mat_tgt,vec_tgt,regparam);
  }

  #   4. preliminary rank determination
  diagN = diag(n)
  M     = t(diagN-W)%*%(diagN-W)
  St    = (t(pX)%*%M%*%pX) + (t(pX)%*%pX)
  r     = as.integer(Matrix::rankMatrix(St))
  if (r < ndim){
    message("* do.cnpe : intrinsic rank of matrix St is smaller than 'ndim'.")
    ndim = r
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION FOR CNPE
  #   1. EVD for t(Xtilde)%*%Xtilde
  #      select Vr and vecSig1
  Xtilde   = t(pX)%*%cbind(t(diagN-W),diagN) #  (D x 2N)
  Xcost    = t(Xtilde)%*%Xtilde              # (2N x 2N)
  eigXcost = base::eigen(Xcost)

  Vr      = eigXcost$vectors[,1:r]
  vecSig1 = as.vector(eigXcost$values[1:r])

  #   2. compute Ur and Sctilde
  invSig1half = diag(1/sqrt(vecSig1))
  Ur          = Xtilde%*%Vr%*%invSig1half
  Sctilde     = invSig1half%*%t(Ur)%*%t(pX)%*%pX%*%Ur%*%invSig1half

  #   3. decompose Sctilde and denote it as Wmat
  Wmat = base::eigen(Sctilde)$vectors

  #   4. use first ndim unitary vectors
  resmat     = Ur%*%invSig1half%*%Wmat
  projection = aux.adjprojection(resmat[,1:ndim])

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
