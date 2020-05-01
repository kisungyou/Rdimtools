#' Kernel Principal Component Analysis
#'
#' Kernel principal component analysis (KPCA/Kernel PCA) is a nonlinear extension of classical
#' PCA using techniques called \href{https://en.wikipedia.org/wiki/Kernel_method}{kernel trick},
#' a common method of introducing nonlinearity by transforming, usually, covariance structure or
#' other gram-type estimate to make it flexible in Reproducing Kernel Hilbert Space.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param kernel a vector containing name of a kernel and corresponding parameters. See also \code{\link{aux.kernelcov}} for complete description of Kernel Trick.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{vars}{variances of projected data / eigenvalues from kernelized covariance matrix.}
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
#' ## try out different settings
#' output1 <- do.kpca(X)                         # default setting
#' output2 <- do.kpca(X,kernel=c("gaussian",5))  # gaussian kernel with large bandwidth
#' output3 <- do.kpca(X,kernel=c("laplacian",1)) # laplacian kernel
#'
#' ## visualize three different projections
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(output1$Y, col=label, pch=19, main="Gaussian kernel")
#' plot(output2$Y, col=label, pch=19, main="Gaussian kernel with sigma=5")
#' plot(output3$Y, col=label, pch=19, main="Laplacian kernel")
#' par(opar)
#'
#' @seealso \code{\link{aux.kernelcov}}
#' @references
#' \insertRef{scholkopf_kernel_1997}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_KPCA
#' @concept nonlinear_methods
#' @export
do.kpca <- function(X,ndim=2,preprocess=c("null","center","scale","cscale","whiten","decorrelate"),kernel=c("gaussian",1.0)){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.kpca : 'ndim' is a positive integer in [1,#(covariates)].")
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
  Koriginal = Ks$K
  Kcentered = Ks$Kcenter

  # 4. main computation
  output = aux.eigendec(Kcentered)
  eigvals = output$eigval
  eigvecs = output$eigvec

  # 5. result
  result = list()
  result$Y = t(t(eigvecs[,1:k]) %*% Kcentered)
  result$trfinfo = trfinfo
  result$vars = (eigvals[1:k])/N
  return(result)
}
