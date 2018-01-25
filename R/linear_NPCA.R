#' Nonnegative Principal Component Analysis
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center" and other options of "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
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
#' @references
#' \insertRef{zafeiriou_nonnegative_2010}{Rdimtools}
#'
#' @rdname linear_NPCA
#' @author Kisung You
#' @export
do.npca <- function(X, ndim=2, preprocess=c("center","decorrelate","whiten"),
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
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR NONNEGATIVE PCA
  #   1. initialize for U
  Uinit = matrix(runif(p*ndim),nrow=p)
  #   2. compute C
  C = cov(pX)
  #   3. compute projection matrix
  projection = aux.adjprojection(method_nnprojmax(C, Uinit, reltol, maxiter))


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
