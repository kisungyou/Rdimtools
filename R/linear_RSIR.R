#' Regularized Sliced Inverse Regression
#'
#' One of possible drawbacks in SIR method is that for high-dimensional data, it might suffer from
#' rank deficiency of scatter/covariance matrix. Instead of naive matrix inversion, several have
#' proposed regularization schemes that reflect several ideas from various incumbent methods.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param response a length-\eqn{n} vector of response variable.
#' @param ndim an integer-valued target dimension.
#' @param h the number of slices to divide the range of response vector.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center" and other options "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param regmethod type of regularization scheme to be used.
#' @param tau regularization parameter for adjusting rank-deficient scatter matrix.
#' @param numpc number of principal components to be used in intermediate dimension reduction scheme.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
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
#' ## try with different regularization methods
#' ## use default number of slices
#' out1 = do.rsir(X, y, regmethod="Ridge")
#' out2 = do.rsir(X, y, regmethod="Tikhonov")
#' out3 = do.rsir(X, y, regmethod="PCA")
#' out4 = do.rsir(X, y, regmethod="PCARidge")
#' out5 = do.rsir(X, y, regmethod="PCATikhonov")
#' outsir = do.sir(X, y)
#'
#' ## visualize
#' par(mfrow=c(2,3))
#' plot(out1$Y[,1], out1$Y[,2], main="RSIR::Ridge")
#' plot(out2$Y[,1], out2$Y[,2], main="RSIR::Tikhonov")
#' plot(out3$Y[,1], out3$Y[,2], main="RSIR::PCA")
#' plot(out4$Y[,1], out4$Y[,2], main="RSIR::PCA+Ridge")
#' plot(out5$Y[,1], out5$Y[,2], main="RSIR::PCA+Tikhonov")
#' plot(outsir$Y[,1], outsir$Y[,2], main="standard SIR")
#'
#' @references
#' \insertRef{chiaromonte_dimension_2002}{Rdimtools}
#'
#' \insertRef{zhong_rsir:_2005}{Rdimtools}
#'
#' \insertRef{bernard-michel_gaussian_2009}{Rdimtools}
#'
#' \insertRef{bernard-michel_retrieval_2009}{Rdimtools}
#'
#'
#' @seealso \code{\link{do.sir}}
#' @author Kisung You
#' @rdname linear_RSIR
#' @export
do.rsir <- function(X, response, ndim=2, h=max(2, round(nrow(X)/5)), preprocess=c("center","decorrelate","whiten"),
                    regmethod=c("Ridge","Tikhonov","PCA","PCARidge","PCATikhonov"), tau=1.0, numpc=ndim){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. response
  response = as.double(response)
  if ((any(is.infinite(response)))||(!is.vector(response))||(any(is.na(response)))){
    stop("* do.rsir : 'response' should be a vector containing no NA values.")
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
  #   6. regmethod
  if (missing(regmethod)){
    regalgorithm = "Ridge"
  } else {
    regalgorithm = match.arg(regmethod)
  }
  #   7. tau
  tau = as.double(tau)
  if (!check_NumMM(tau,0,Inf,compact=FALSE)){stop("* do.rsir : 'tau' should be a strictly positive real number.")}
  #   8. numpc
  numpc = as.integer(numpc)

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
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
  #   1. take the case branching
  if (regalgorithm=="Ridge"){
    Omega = (1/tau)*diag(p)
  } else if (regalgorithm=="Tikhonov"){
    Omega = (1/tau)*mat_Sigma
  } else {
    pcarun = RSpectra::eigs(mat_Sigma,numpc)
    pcavalues = pcarun$values
    pcavector = pcarun$vectors
    if (regalgorithm=="PCA"){
      Omega = array(0,c(p,p))
      for (i in 1:numpc){
        cvec  = as.vector(pcavector[,i])
        Omega = Omega + (1/as.double(pcavalues[i]))*outer(cvec,cvec)
      }
    } else if (regalgorithm=="PCARidge"){
      Omega = array(0,c(p,p))
      for (i in 1:numpc){
        cvec  = as.vector(pcavector[,i])
        Omega = Omega + (1/tau)*outer(cvec,cvec)
      }
    } else if (regalgorithm=="PCATikhonov"){
      Omega = array(0,c(p,p))
      for (i in 1:numpc){
        cvec  = as.vector(pcavector[,i])
        Omega = Omega + (1/tau)*(as.double(pcavalues[i]))*outer(cvec,cvec)
      }
    }
  }
  #   2. adjust LHS and RHS terms
  LHS = (Omega%*%mat_Sigma)+diag(p)
  RHS = (Omega%*%mat_Gamma)
  #   2. perform matrix inversion
  costInv = aux.bicgstab(LHS, RHS, verbose=FALSE)$x
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
