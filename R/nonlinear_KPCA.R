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
#' Default is ``null'', and other methods of ``decorrelate'',``center'' , and ``whiten'' are supported. See also \code{\link{aux.preprocess}} for more details.
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
#' \dontrun{
#' ## generate ribbon-shaped data
#' X = aux.gensamples(dname="ribbon",n=123)
#'
#' ## 1. standard KPCA with gaussian kernel
#' output1 <- do.kpca(X,ndim=2)
#'
#' ## 2. gaussian kernel with large bandwidth
#' output2 <- do.kpca(X,ndim=2,kernel=c("gaussian",5))
#'
#' ## 3. use laplacian kernel
#' output3 <- do.kpca(X,ndim=2,kernel=c("laplacian",1))
#'
#' ## Visualize three different projections
#' par(mfrow=c(1,3))
#' plot(output1$Y[,1],output1$Y[,2],main="Gaussian kernel")
#' plot(output2$Y[,1],output2$Y[,2],main="Gaussian kernel with sigma=5")
#' plot(output3$Y[,1],output3$Y[,2],main="Laplacian kernel")
#' }
#'
#' @seealso \code{\link{aux.kernelcov}}
#' @references
#' \insertRef{goos_kernel_1997}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_KPCA
#' @export
do.kpca <- function(X,ndim=2,preprocess="null",kernel=c("gaussian",1.0)){
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
  algpreprocess = preprocess
  if (!is.element(algpreprocess,c("null","center","whiten","decorrelate"))){
    stop("* do.kpca : 'preprocess' argument is invalid.")
  }
  ktype = kernel

  # 3. preprocess
  #   3-1. centering or so.
  if (preprocess=="null"){
    trfinfo = list()
    trfinfo$type = "null"
    pX = as.matrix(X,nrow=nrow(X))
  } else {
    tmplist = aux.preprocess(X,type=preprocess)
    trfinfo = tmplist$info
    pX      = tmplist$pX
  }
  trfinfo$algtype = "nonlinear"

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
