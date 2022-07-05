#' Forward Orthogonal Search by Maximizing the Overall Dependency
#'
#' The FOS-MOD algorithm \insertCite{wei_2007_FeatureSubsetSelection}{Rdimtools}
#' is an unsupervised algorithm that selects a desired number of features in
#' a forward manner by ranking the features using the squared correlation
#' coefficient and sequential orthogonalization.
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
#' out1 = do.fosmod(iris.dat)
#' out2 = do.lscore(iris.dat)
#' out3 = do.fscore(iris.dat, iris.lab)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=iris.lab, main="FOS-MOD")
#' plot(out2$Y, pch=19, col=iris.lab, main="Laplacian Score")
#' plot(out3$Y, pch=19, col=iris.lab, main="Fisher Score")
#' par(opar)
#' }
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param ... extra parameters including \describe{
#' \item{preprocess}{an additional option for preprocessing the data.
#' See also \code{\link{aux.preprocess}} for more details (default: \code{"center"}).}
#' }
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @rdname feature_FOSMOD
#' @concept feature_methods
#' @export
do.fosmod <- function(X, ndim=2, ...){
  #------------------------------------------------------------------------
  # INPUT : EXPLICIT
  # data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)

  # ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.fosmod : 'ndim' is a positive integer in [1,#(covariates)].")
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

  #------------------------------------------------------------------------
  # COMPUTATION : PRELIMINARY
  tmplist = aux.preprocess.hidden(X, type=par_preprocess, algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  # COMPUTATION : INITIAL
  Cbar1     = as.vector(base::rowMeans(base::abs(stats::cor(pX))))
  now_id    = which.max(Cbar1)[1]
  now_ortho = as.vector(pX[,now_id])

  #------------------------------------------------------------------------
  # COMPUTATION : ITERATIVE
  for (it in 2:ndim){
    # target IDs
    target_id = setdiff(1:p, now_id)

    # target vectors orthogonalized
    if (it < 3){
      target_Qs = cpp_fosmod_orthogonalize_vec(now_ortho, pX[,target_id])
    } else {
      target_Qs = cpp_fosmod_orthogonalize(now_ortho, pX[,target_id])
    }

    # compute the squared-correlation coefficient vector
    target_cost = as.vector(cpp_fosmod_crosscorr(pX, target_Qs))

    # find the best (in terms of the available IDs)
    best_id = which.max(target_cost)

    # update the information
    now_id    = c(now_id, target_id[best_id])
    if (it < 3){
      now_ortho = cbind(matrix(now_ortho, ncol=1), matrix(as.vector(target_Qs[,best_id]), ncol=1))
    } else {
      now_ortho = cbind(now_ortho, as.vector(target_Qs[,best_id]))
    }
  }

  # find the projection
  projection = aux.featureindicator(p,ndim,now_id)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$featidx = now_id
  result$projection = projection
  result$trfinfo    = trfinfo
  result$algorithm  = "linear:FOSMOD"
  return(structure(result, class="Rdimtools"))
}



