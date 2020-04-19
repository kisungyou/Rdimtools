#' Supervised Principal Component Analysis
#'
#' Unlike original principal component analysis (\code{\link{do.pca}}), this algorithm implements
#' a supervised version using response information for feature selection. For each feature/column,
#' its normalized association with \code{response} variable is computed and the features with
#' large magnitude beyond \code{threshold} are selected. From the selected submatrix,
#' regular PCA is applied for dimension reduction.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param response a length-\eqn{n} vector of response variable.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is \code{center}. See also \code{\link{aux.preprocess}} for more details.
#' @param threshold a threshold value to cut off normalized association between covariates and response.
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
#' ## try different threshold values
#' out1 = do.spc(X, y, threshold=2)
#' out2 = do.spc(X, y, threshold=5)
#' out3 = do.spc(X, y, threshold=10)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="SPC::threshold=2")
#' plot(out2$Y, main="SPC::threshold=5")
#' plot(out3$Y, main="SPC::threshold=10")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{bair_prediction_2006}{Rdimtools}
#'
#' @rdname linear_SPC
#' @author Kisung You
#' @concept linear_methods
#' @export
do.spc <- function(X, response, ndim=2, preprocess=c("center","whiten","decorrelate"), threshold=0.1){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. response
  response = as.numeric(response)
  if ((any(is.infinite(response)))||(!is.vector(response))||(any(is.na(response)))||(length(response)!=n)){
    stop("* do.spc : 'response' should be a vector containing no NA values.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.spc : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. threshold
  if ((length(threshold)>1)||(threshold<=0)){
    stop("* do.spc : 'threshold' must be a positive real number.")
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : DATA PREPROCESSING
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## Main Computation
  #   1. find the ones with bigger univariate measurements
  svec = rep(0,p)
  for (j in 1:p){
    xj      = as.vector(pX[,j])
    svec[j] = sum(xj*response)/sqrt(sum(xj*xj))
  }
  uidx = which(abs(svec)>threshold)
  if (length(uidx)<ndim){
    ordsvec = sort(abs(svec), decreasing=FALSE)
    uidx    = which(abs(svec)>=ordsvec[ndim])
  }
  pXpart    = pX[,uidx]
  #   2. compute largest eigenvectors
  projection = aux.adjprojection(RSpectra::eigs(stats::cov(pXpart), ndim, which="LM")$vectors)

  #------------------------------------------------------------------------
  ## Return Results
  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = (pXpart%*%projection)
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}

