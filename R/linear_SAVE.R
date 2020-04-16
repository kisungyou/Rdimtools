#' Sliced Average Variance Estimation
#'
#' Sliced Average Variance Estimation (SAVE) is a supervised linear dimension reduction method.
#' It is based on sufficiency principle with respect to central subspace concept under the
#' linerity and constant covariance conditions. For more details, see the reference paper.
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
#' out1 = do.save(X, y, h=2)
#' out2 = do.save(X, y, h=5)
#' out3 = do.save(X, y, h=10)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="SAVE::2 slices")
#' plot(out2$Y, main="SAVE::5 slices")
#' plot(out3$Y, main="SAVE::10 slices")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{denniscook_method_2000}{Rdimtools}
#'
#' @seealso \code{\link{do.sir}}
#' @author Kisung You
#' @rdname linear_SAVE
#' @concept linear_methods 
#' @export
do.save <- function(X, response, ndim=2, h=max(2, round(nrow(X)/5)),
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
  if (!check_ndim(ndim,p)){stop("* do.save : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. h : number of slices
  h = as.integer(h)
  if (!is.factor(response)){
    if (!check_NumMM(h,2,ceiling(n/2),compact=TRUE)){stop("* do.save : the number of slices should be in [2,n/2].")}
  }
  #   5. preprocess
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
  #   3. compute scaling as noted from the method.
  #   3-1. mean
  tmp_mean = colMeans(pX)
  #   3-2. eigendecomposition of covariance matrix
  tmp_SigmaInvHalf = save_SigmaInvHalf(pX)
  #   3-3. adjust values
  tmpZ = array(0,c(n,p))
  for (i in 1:n){
    tmpZ[i,] = as.vector(pX[i,])-tmp_mean
  }
  Z    = tmpZ%*%tmp_SigmaInvHalf
  #   4. Construct M
  M = array(0,c(p,p))
  for (s in 1:nlabel){
    idxs = which(label==ulabel[s])
    ns   = length(idxs)
    IVs  = diag(p)-cov(pX[idxs,]) # difference !
    M    = M + (ns/n)*(IVs%*%IVs)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION
  projection = aux.adjprojection(RSpectra::eigs(M, ndim)$vectors)

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
save_SigmaInvHalf <- function(X){
  coveig = eigen(cov(X))
  # 1. adjust eigenvalues
  evalues = coveig$values
  nevals  = length(evalues)
  newvals = rep(0,nevals)
  for (i in 1:nevals){
    tgt = evalues[i]
    if (tgt>0){
      newvals[i] = 1/sqrt(tgt)
    }
  }
  # 2. eigenvectors
  Mlambda = coveig$vectors
  # 3. compute output
  output = Mlambda%*%diag(newvals)%*%t(Mlambda)
  return(output)
}
