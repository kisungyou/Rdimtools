#' Sliced Inverse Regression
#'
#' Sliced Inverse Regression (SIR) is a supervised linear dimension reduction technique.
#' Unlike engineering-driven methods, SIR takes a concept of \emph{central subspace}, where
#' conditional independence after projection is guaranteed. It first divides the range of
#' response variable. Projection vectors are extracted where projected data best explains response variable.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param response a length-\eqn{n} vector of response variable.
#' @param ndim an integer-valued target dimension.
#' @param h the number of slices to divide the range of response vector.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \donttest{
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
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="SIR::2 slices")
#' plot(out2$Y, main="SIR::5 slices")
#' plot(out3$Y, main="SIR::10 slices")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{li_sliced_1991}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_SIR
#' @concept linear_methods 
#' @export
do.sir <- function(X, response, ndim=2, h=max(2, round(nrow(X)/5)),
                   preprocess=c("center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. response
  response = as.double(response)
  if ((any(is.infinite(response)))||(!is.vector(response))||(any(is.na(response)))){
    stop("* do.sir : 'response' should be a vector containing no NA values.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.sir : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. h : number of slices
  h = as.integer(h)
  if (!is.factor(response)){
    if (!check_NumMM(h,2,ceiling(n/2),compact=TRUE)){stop("* do.save : the number of slices should be in [2,n/2].")}
  }  #   5. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #   2. build label matrix
  if (!is.factor(response)){
    label  = as.integer(sir_makelabel(response, h))
  } else {
    label = as.integer(response)
  }
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
  costInv = aux.bicgstab(mat_Sigma, mat_Gamma, verbose=FALSE)$x
  #   2. find top eigenvectors
  projection = aux.adjprojection(RSpectra::eigs(costInv, ndim)$vectors)
  projection = matrix(as.double(projection), nrow=p)


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

