#' Kernel Entropy Component Analysis
#'
#' Kernel Entropy Component Analysis(KECA) is a kernel method of dimensionality reduction.
#' Unlike Kernel PCA(\code{\link{do.kpca}}), it utilizes eigenbasis of kernel matrix \eqn{K}
#' in accordance with indices of largest Renyi quadratic entropy in which entropy for
#' \eqn{j}-th eigenpair is defined to be \eqn{\sqrt{\lambda_j}e_j^T 1_n}, where \eqn{e_j} is
#' \eqn{j}-th eigenvector of an uncentered kernel matrix \eqn{K}.
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param kernel a vector containing name of a kernel and corresponding parameters. See also \code{\link{aux.kernelcov}} for complete description of Kernel Trick.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{entropy}{a length-\code{ndim} vector of estimated entropy values.}
#' }
#'
#' @examples
#' ## load iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' ## 1. standard KECA with gaussian kernel
#' output1 <- do.keca(X,ndim=2)
#'
#' ## 2. gaussian kernel with large bandwidth
#' output2 <- do.keca(X,ndim=2,kernel=c("gaussian",5))
#'
#' ## 3. use laplacian kernel
#' output3 <- do.keca(X,ndim=2,kernel=c("laplacian",1))
#'
#' ## Visualize three different projections
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(output1$Y, pch=19, col=label, main="Gaussian kernel")
#' plot(output2$Y, pch=19, col=label, main="Gaussian, sigma=5")
#' plot(output3$Y, pch=19, col=label, main="Laplacian kernel")
#' par(opar)
#'
#' @seealso \code{\link{aux.kernelcov}}
#' @references
#' \insertRef{jenssen_kernel_2010}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_KECA
#' @concept nonlinear_methods
#' @export
do.keca <- function(X,ndim=2,kernel=c("gaussian",1.0),
                    preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.keca : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  k = as.integer(ndim)
  N = nrow(X)
  d = ncol(X)
  # 2. ... parameters
  #   preprocess : 'null', 'center','decorrelate', or 'whiten'
  #   kernel     : c("kernel name",par1,par2)
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  ktype = kernel

  # 3. preprocess
  #   3-1. centering or so.
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   3-2. compute K and centered K
  Ks = aux.kernelcov(pX,ktype)
  K  = Ks$K

  # 4. main computation
  #   4-1 eigendecomposition
  outeigen = eigen(K)
  evals    = outeigen$values
  evecs    = outeigen$vectors
  onesn    = as.matrix(array(1,c(nrow(evecs),1))) # ones(N,1)

  #   4-2. entropy computation
  entvec = array(0,c(1,length(evals)))
  for (i in 1:length(entvec)){
    entvec[i] = ifelse(evals[i]>=0,sqrt(evals[i])*sum(evecs[,i]*onesn),-10000)
  }
  largestidx = order(entvec,decreasing=TRUE)[1:ndim]

  #   4-2. results
  outputentvec = entvec[largestidx]
  Y = t(diag(sqrt(evals[largestidx]))%*%t(evecs[,largestidx]))
  result = list()
  result$Y = Y
  result$trfinfo = trfinfo
  result$entropy = outputentvec
  return(result)
}
