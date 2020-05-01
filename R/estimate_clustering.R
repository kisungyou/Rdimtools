#' Intrinsic Dimension Estimation via Clustering
#'
#' Instead of directly using neighborhood information, \code{est.clustering} adopts hierarchical
#' neighborhood information using \code{\link[stats]{hclust}} by recursively merging leafs
#' over the range of radii.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param kmin minimal number of neighborhood size to search over.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{estimated intrinsic dimension.}
#' }
#'
#' @examples
#' \donttest{
#' ## create 'swiss' roll dataset
#' X = aux.gensamples(dname="swiss")
#'
#' ## try different k values
#' out1 = est.clustering(X, kmin=5)
#' out2 = est.clustering(X, kmin=25)
#' out3 = est.clustering(X, kmin=50)
#'
#' ## print the results
#' line1 = paste0("* est.clustering : kmin=5  gives ",round(out1$estdim,2))
#' line2 = paste0("* est.clustering : kmin=25 gives ",round(out2$estdim,2))
#' line3 = paste0("* est.clustering : kmin=50 gives ",round(out3$estdim,2))
#' cat(paste0(line1,"\n",line2,"\n",line3))
#' }
#'
#' @references
#' \insertRef{eriksson_estimating_2012}{Rdimtools}
#'
#' @rdname estimate_clustering
#' @author Kisung You
#' @export
est.clustering <- function(X, kmin=round(sqrt(nrow(X)))){
  ##########################################################################
  ## Preprocessing and Default Parameter
  aux.typecheck(X)
  n  = nrow(X)
  D  = as.matrix(dist(X))

  rmax = max(D)
  rmin = min(apply(D,1,sort)[kmin+1,])

  ##########################################################################
  ## Main Body
  #   1. construct That using complete agglomerative clustering
  That     = stats::hclust(as.dist(D), method="complete")
  clusters = stats::cutree(That, h=That$height)
  histcl   = array(0,dim(clusters))
  #   2. iteration
  for (i in 1:ncol(clusters)){
    cclust  = clusters[,i]
    for (j in 1:length(unique(cclust))){
      idxj = (cclust==j)
      histcl[idxj,i] = max(D[idxj,idxj])
    }
  }
  #   3. ID estimation via lm fitting
  nslice = round(n/4)
  vecr   = seq(from=rmin,to=rmax,length.out=nslice)
  logr   = log(vecr)
  logM   = rep(0,nslice)
  for (i in 1:nslice){
    logM[i] = log(length(unique(histcl[, which(colMeans(histcl > vecr[i]) == 1)[1]])))
  }
  Dpack = abs(as.numeric(coefficients(lm(logM~logr))[2]))


  ##########################################################################
  ## Result plot
  result = list()
  result$estdim = max(Dpack,1)
  return(result)
}
