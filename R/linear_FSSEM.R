#' Feature Subset Selection using Expectation-Maximization
#'
#' Feature Subset Selection using Expectation-Maximization (FSSEM) takes a wrapper approach to feature selection problem.
#' It iterates over optimizing the selection of variables by incrementally including each variable that adds the most
#' significant amount of scatter separability from a labeling obtained by Gaussian mixture model. This method is
#' quite computation intensive as it pertains to multiple fitting of GMM. Setting smaller \code{max.k} for each round of
#' EM algorithm as well as target dimension \code{ndim} would ease the burden.
#'
#' @examples
#' ## run FSSEM with IRIS dataset - select 2 of 4 variables
#' data(iris)
#' irismat = as.matrix(iris[,2:4])
#'
#' ## select 50 observations for CRAN-purpose small example
#' id50 = sample(1:nrow(irismat), 50)
#' sel.dat = irismat[id50,]
#' sel.lab = as.factor(iris[id50,5])
#'
#' ## run and visualize
#' out0 = do.fssem(sel.dat, ndim=2, max.k=3)
#' opar = par(no.readonly=TRUE)
#' plot(out0$Y, main="small run", col=sel.lab, pch=19)
#' par(opar)
#'
#' \dontrun{
#' ## NOT-FOR-CRAN example; run at your machine !
#' ## try different maximum number of clusters
#' out3 = do.fssem(irismat, ndim=2, max.k=3)
#' out6 = do.fssem(irismat, ndim=2, max.k=6)
#' out9 = do.fssem(irismat, ndim=2, max.k=9)
#'
#' ## visualize
#' cols = as.factor(iris[,5])
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(3,1))
#' plot(out3$Y, main="max k=3", col=cols)
#' plot(out6$Y, main="max k=6", col=cols)
#' plot(out9$Y, main="max k=9", col=cols)
#' par(opar)
#' }
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param max.k maximum number of clusters for GMM fitting with EM algorithms.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @references
#' \insertRef{dy_feature_2004}{Rdimtools}
#'
#' @rdname linear_FSSEM
#' @author Kisung You
#' @export
do.fssem <- function(X, ndim=2, max.k=10, preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.fssem : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. max.k : maximum number of potential clusters
  max.k = round(max.k)
  if (max.k > round(nrow(X)/2)){
    stop("* do.fssem : 'max.k' should be a reasonable upper bound for potential clusters number.")
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #  preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION : ITERATE UNTIL NDIM
  idalls  = seq(from=1,to=p,by=1)
  idgood  = c()
  for (iter in 1:ndim){
    # 1. generate a list
    if (iter > 1){
      addnow = base::setdiff(idalls, idgood)
    } else {
      addnow = idalls
    }
    nnnnow = length(addnow)

    # 2. prepare to record
    scores = rep(0, nnnnow)

    # 3. do the iteration, now !
    for (i in 1:nnnnow){
      # 3-1. current set of ids
      currentids = c(idgood, addnow[i])
      # 3-2. compute the EM-based partition
      cpartition = fssem.single(pX[,currentids], max.k)
      # 3-3. compute scatter separability
      if (length(unique(cpartition))<2){ # it's possible that 1-partition is the best.. ignore !
        scores[i] = 0
      } else {
        scatmat = aux.2scatter(pX[,currentids], cpartition)
        if (iter<2){
          scores[i] = scatmat$between/scatmat$within
        } else {
          scores[i] = sum(diag(aux.pinv(scatmat$within)%*%(scatmat$between)))
        }
      }
    }

    # 4. select the best
    wmscore = which.max(scores) # id of the best !
    if (length(wmscore)>1){
      wmscore = wmscore[1]
    }
    idgood = c(idgood, addnow[wmscore])
  }
  #   4. find the projection matrix
  projection = aux.featureindicator(p,ndim,idgood)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}



#   -----------------------------------------------------------------------
#' @keywords internal
#' @noRd
fssem.single <- function(Xpart, maxk){ # return the optimal value
  if (is.vector(Xpart)){
    Xpart = matrix(Xpart)
  }
  xdist  = stats::dist(Xpart)
  myfun  = utils::getFromNamespace("hidden_kmedoids_best","maotai")
  hey    = myfun(xdist, mink=1, maxk=round(maxk))
  return(as.vector(hey$label[,hey$opt.k]))

  # # 1. fit GMM
  # bicvalues = ClusterR::Optimal_Clusters_GMM(Xpart, maxk, criterion="BIC", verbose=FALSE, plot_data=FALSE)
  # # 2. select out the best
  # koptimal  = which.min(bicvalues)
  # # 3. run GMM
  # emoutput = ClusterR::GMM(Xpart, gaussian_comps = koptimal)
  # # 4. cluster label
  # pr = as.vector(ClusterR::predict_GMM(Xpart, emoutput$centroids, emoutput$covariance_matrices, emoutput$weights)$cluster_labels)
  # return(pr)
}
