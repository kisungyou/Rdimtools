#' Nonnegative Principal Component Analysis
#'
#' Nonnegative Principal Component Analysis (NPCA) is a variant of PCA where
#' projection vectors - or, basis for learned subspace - contain no negative values.
#'
#' @param X an \eqn{(n\times p)} matrix whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations (default: 100).}
#' \item{reltol}{relative tolerance stopping criterion (default: 1e-4).}
#' }
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are principal components.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' \dontrun{
#' ## use iris data
#' data(iris, package="Rdimtools")
#' set.seed(100)
#' subid = sample(1:150, 50)
#' X     = as.matrix(iris[subid,1:4]) + 50
#' label = as.factor(iris[subid,5])
#'
#' ## run NCPA and compare with others
#' outNPC = do.npca(X)
#' outPCA = do.pca(X)
#' outMVP = do.mvp(X, label)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(outNPC$Y, pch=19, col=label, main="NPCA")
#' plot(outPCA$Y, pch=19, col=label, main="PCA")
#' plot(outMVP$Y, pch=19, col=label, main="MVP")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zafeiriou_nonnegative_2010}{Rdimtools}
#'
#' @seealso \code{\link{do.pca}}
#' @rdname linear_NPCA
#' @author Kisung You
#' @concept linear_methods
#' @export
do.npca <- function(X, ndim=2, ...){
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
  # if (missing(preprocess)){
  #   algpreprocess = "center"
  # } else {
  #   algpreprocess = match.arg(preprocess)
  # }
  #   * maxiter and reltol

  # Extra parameters
  params  = list(...)
  pnames  = names(params)

  if ("reltol"%in%pnames){
    reltol = max(.Machine$double.eps, as.double(params$reltol))
  } else {
    reltol = 10^(-4)
  }
  if ("maxiter"%in%pnames){
    maxiter = max(5, round(params$maxiter))
  } else {
    maxiter = 100
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : DATA PREPROCESSING
  # tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  # trfinfo = tmplist$info
  # pX      = tmplist$pX
  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR NONNEGATIVE PCA
  #   1. initialize for U
  Uinit = matrix(runif(p*ndim),nrow=p)
  #   2. compute C
  C = stats::cov(X)
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
  result$Y = X%*%projection
  result$projection = projection
  result$algorithm  = "linear:NPCA"
  return(structure(result, class="Rdimtools"))
}
