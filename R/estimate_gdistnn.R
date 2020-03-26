#' Intrinsic Dimension Estimation based on Manifold Assumption and Graph Distance
#'
#' As the name suggests, this function assumes that the data is sampled from the manifold in that
#' graph representing the underlying manifold is first estimated via \eqn{k}-nn. Then graph distance
#' is employed as an approximation of geodesic distance to locally estimate intrinsic dimension.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param k the neighborhood size used for constructing a graph. We suggest it to be large enough to build a connected graph.
#' @param k1 local neighborhood parameter (smaller radius) for graph distance.
#' @param k2 local neighborhood parameter (larger radius)  for graph distance.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{the global estimated dimension, which is averaged local dimension.}
#' \item{estloc}{a length-\eqn{n} vector of locally estimated dimension at each point.}
#' }
#'
#' @examples
#' \donttest{
#' ## create 3 datasets of intrinsic dimension 2.
#' X1 = aux.gensamples(dname="swiss")
#' X2 = aux.gensamples(dname="ribbon")
#' X3 = aux.gensamples(dname="saddle")
#'
#' ## acquire an estimate for intrinsic dimension
#' out1 = est.gdistnn(X1, k=10)
#' out2 = est.gdistnn(X2, k=10)
#' out3 = est.gdistnn(X3, k=10)
#'
#' ## print the results
#' sprintf("* est.gdistnn : estimated dimension for 'swiss'  data is %.2f.",out1$estdim)
#' sprintf("* est.gdistnn : estimated dimension for 'ribbon' data is %.2f.",out2$estdim)
#' sprintf("* est.gdistnn : estimated dimension for 'saddle' data is %.2f.",out3$estdim)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' hist(out1$estloc, main="Result-'Swiss'", xlab="local dimension")
#' abline(v=out1$estdim, lwd=3, col="red")
#' hist(out2$estloc, main="Result-'Ribbon'", xlab="local dimension")
#' abline(v=out2$estdim, lwd=3, col="red")
#' hist(out3$estloc, main="Result-'Saddle'", xlab="local dimension")
#' abline(v=out2$estdim, lwd=3, col="red")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{he_intrinsic_2014}{Rdimtools}
#'
#' @rdname estimate_gdistnn
#' @author Kisung You
#' @export
est.gdistnn <- function(X, k=5, k1=3, k2=10){
  ##########################################################################
  ## preprocessing
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)

  myk  = round(k)
  myk1 = round(k1)
  myk2 = round(k2)

  ##########################################################################
  ## computation 1 : we want a connected graph
  nbdtype   = c("knn",myk)
  nbdstruct = aux.graphnbd(X,method="euclidean",
                           type=nbdtype,symmetric="union")
  D     = nbdstruct$dist
  Dmask = nbdstruct$mask
  nD    = ncol(D)
  wD    = Dmask*D
  wD[is.na(wD)] = 0.0
  sD = aux.shortestpath(wD)
  if (any(is.infinite(sD))){
    stop("* est.gdistnn : a graph is not connected. Instead of extracting maximal subgraph, we propose to increase the value of 'k'.")
  }

  ##########################################################################
  ## computation 2 : per-sample computation
  #   2-1. compute k1-th and k2-th values
  vals.k1 = rep(0,n)
  vals.k2 = rep(0,n)
  for (i in 1:n){
    tgtvec = base::sort(as.vector(sD[,i]))
    vals.k1[i] = tgtvec[1+myk1]
    vals.k2[i] = tgtvec[1+myk2]
  }
  #   2-2. local values
  vec.output = rep(0,n)
  for (i in 1:n){
    term1 = log(k1)-log(k2)
    term2 = log(vals.k1[i]) - log(vals.k2[i])
    vec.output[i] = term1/term2
  }
  #   2-3. global value
  val.output = base::mean(vec.output)


  ##########################################################################
  ## Return the results
  result = list()
  result$estdim = val.output
  result$estloc = vec.output
  return(result)
}
