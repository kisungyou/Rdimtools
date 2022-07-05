#' Principal Feature Analysis
#'
#' Principal Feature Analysis \insertCite{lu_2007_FeatureSelectionUsing}{Rdimtools}
#' adopts an idea from the celebrated PCA for unsupervised feature selection.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param ... extra parameters including \describe{
#' \item{cor}{mode of eigendecomposition. \code{FALSE} for decomposing the
#' empirical covariance matrix and \code{TRUE} uses the correlation matrix
#' (default: \code{FALSE}).}
#' \item{preprocess}{an additional option for preprocessing the data.
#' See also \code{\link{aux.preprocess}} for more details (default: \code{"center"}).}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' ## it is known that feature 3 and 4 are more important.
#' data(iris)
#' set.seed(100)
#' subid    <- sample(1:150, 50)
#' iris.dat <- as.matrix(iris[subid,1:4])
#' iris.lab <- as.factor(iris[subid,5])
#'
#' ## compare with other methods
#' out1 = do.pfa(iris.dat)
#' out2 = do.lscore(iris.dat)
#' out3 = do.fscore(iris.dat, iris.lab)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=iris.lab, main="Principal Feature Analysis")
#' plot(out2$Y, pch=19, col=iris.lab, main="Laplacian Score")
#' plot(out3$Y, pch=19, col=iris.lab, main="Fisher Score")
#' par(opar)
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @rdname feature_PFA
#' @concept feature_methods
#' @export
do.pfa <- function(X, ndim=2, ...){
  #------------------------------------------------------------------------
  # INPUT : EXPLICIT
  # data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)

  # ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.pfa : 'ndim' is a positive integer in [1,#(covariates)].")
  }

  # INPUT : IMPLICIT
  params = list(...)
  pnames = names(params)

  # preprocessing
  if ("preprocess"%in%pnames){
    par_preprocess = tolower(params$preprocess)
  } else {
    par_preprocess = "center"
  }
  if (par_preprocess%in%c("null","scale")){
    stop("* do.pfa : PFA does not allow a preprocessing scheme without centering.")
  }
  if ("cor"%in%pnames){
    par_usecor = as.logical(params$cor)
  } else {
    par_usecor = FALSE
  }


  #------------------------------------------------------------------------
  # COMPUTATION : PRELIMINARY
  tmplist = aux.preprocess.hidden(X, type=par_preprocess, algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  # COMPUTATION : MAIN
  # empirical 2nd order matrix
  if (par_usecor){
    S = stats::cor(pX)
  } else {
    S = stats::cov(pX)
  }

  # compute the (ndim-1) largest eigenvectors
  eigS = RSpectra::eigs_sym(S, (ndim-1))
  Aq   = eigS$vectors

  # perform k-means clustering on the row-space
  run_kmeans = stats::kmeans(base::abs(Aq), ndim)

  # compute the cross distance
  cross_dist = v2aux_pdist2(Aq, run_kmeans$centers)

  # for each center, find the nearest
  idxvec = rep(0, ndim)
  for (i in 1:ndim){
    idxvec[i] = which.min(as.vector(cross_dist[,i]))
  }

  # find the projection
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$featidx = idxvec
  result$projection = projection
  result$trfinfo    = trfinfo
  result$algorithm  = "linear:PFA"
  return(structure(result, class="Rdimtools"))
}
