#' Nonnegative Principal Component Analysis
#'
#' Nonnegative Principal Component Analysis (NPCA) is a variant of PCA where
#' projection vectors - or, basis for learned subspace - contain no negative values.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param maxiter number of maximum iteraions allowed.
#' @param reltol stopping criterion for incremental relative error.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])+50
#' label = as.integer(iris$Species)
#'
#' ## use different preprocessing
#' out1 = do.npca(X, preprocess="center")
#' out2 = do.npca(X, preprocess="cscale")
#' out3 = do.npca(X, preprocess="whiten")
#'
#' ## visualize
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(out1$Y, col=label, main="NPCA:: center")
#' plot(out2$Y, col=label, main="NPCA:: cscale")
#' plot(out3$Y, col=label, main="NPCA:: whiten")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zafeiriou_nonnegative_2010}{Rdimtools}
#'
#' @seealso \code{\link{do.pca}}
#' @rdname linear_NPCA
#' @author Kisung You
#' @export
do.npca <- function(X, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"),
                    maxiter=1000, reltol=1e-5){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.npca : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   * maxiter and reltol
  maxiter = as.integer(maxiter)
  if (!check_NumMM(maxiter, 5, 1e+6)){stop("* do.npca : 'maxiter' is a large positive integer for the number of iterations.")}
  reltol = as.double(reltol)
  if (!check_NumMM(reltol, .Machine$double.eps, 1)){stop("* do.npca : 'reltol' is a small positive real number for stopping criterion.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : DATA PREPROCESSING
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR NONNEGATIVE PCA
  #   1. initialize for U
  Uinit = matrix(runif(p*ndim),nrow=p)
  #   2. compute C
  C = cov(pX)
  #   3. compute projection matrix
  projection = method_nnprojmax(C, Uinit, reltol, maxiter)
  #   4. additional step : NA
  projection[(is.na(projection)||(is.infinite(projection)))] = 1
  for (i in 1:ndim){
    tgt = as.vector(projection[,i])
    projection[,i] = tgt/sqrt(sum(tgt^2))
  }
  projection = aux.adjprojection(projection)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
