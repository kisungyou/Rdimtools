#' Sliced Inverse Regression
#'
#'
#'
#' @examples
#' ## generate swiss roll with auxiliary dimensions
#' ## it follows reference example from LSIR paper.
#' n = 123
#' theta = runif(n)
#' h     = runif(n)
#' t     = (1+2*theta)*(3*pi/2)
#' X     = array(0,c(n,10))
#' X[,1] = t*cos(t)
#' X[,2] = 21*h
#' X[,3] = t*sin(t)
#' X[,4:10] = matrix(runif(7*n), nrow=n)
#'
#' ## corresponding response vector
#' y = sin(5*pi*theta)+(runif(n)*sqrt(0.1))
#'
#' ## try with different numbers of slices
#' out1 = do.sir(X, y, h=2)
#' out2 = do.sir(X, y, h=5)
#' out3 = do.sir(X, y, h=10)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="2 slices")
#' plot(out2$Y[,1], out2$Y[,2], main="5 slices")
#' plot(out3$Y[,1], out3$Y[,2], main="10 slices")
#'
#'
#' @author Kisung You
#' @rdname linear_SIR
#' @export
do.sir <- function(X, response, ndim=2, h=max(2, round(nrow(X)/5)), preprocess=c("center","decorrelate","whiten")){
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
  if (!check_ndim(ndim,p)){stop("* do.sir : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. h : number of slices
  h = as.integer(h)
  if (!check_NumMM(h,2,ceiling(n/2),compact=TRUE)){stop("* do.sir : the number of slices should be in [2,n/2].")}
  #   5. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

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
  #   5. compute Between-Slice Covariance
  mat_Gamma = array(0,c(p,p))
  for (i in 1:nlabel){
    vecdiff = (as.vector(class_mean[i,])-(all_mean))
    mat_Gamma = mat_Gamma + outer(vecdiff,vecdiff)*as.double(class_count[i])/n
  }


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
sir_makelabel <- function(responsevec, h){
  n = length(responsevec)
  output = rep(0,n)
  hn = floor(n/h)
  startidx = 1
  currentlabel = 1
  orderlabel = order(responsevec)
  for (i in 1:(h-1)){
    # zeros, find the index
    idxorders = orderlabel[startidx:(startidx+hn-1)]
    # first, fill in the label
    output[idxorders] = currentlabel
    # second, update
    startidx = startidx + hn
    currentlabel = currentlabel + 1
  }
  # final
  idxorders = orderlabel[startidx:n]
  output[idxorders] = currentlabel
  return(output)
}

