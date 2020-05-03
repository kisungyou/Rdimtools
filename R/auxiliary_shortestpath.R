#' Find shortest path using Floyd-Warshall algorithm
#'
#' This is a fast implementation of Floyd-Warshall algorithm to find the
#' shortest path in a pairwise sense using 'RcppArmadillo'. A logical input
#' is also accepted.
#'
#' @param dist either an \eqn{(n\times n)} matrix or a \code{dist} class object.
#' @return an \eqn{(n\times n)} matrix containing pairwise shortest path.
#'
#' @examples
#' \donttest{
#' ## generate a toy data
#' X = aux.gensamples(n=10)
#'
#' ## Find knn graph with k=5
#' Xgraph = aux.graphnbd(X,type=c("knn",5))
#'
#' ## Separately use binarized and real distance matrices
#' W1 = aux.shortestpath(Xgraph$mask)
#' W2 = aux.shortestpath(Xgraph$dist)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' image(W1); title("from binarized")
#' image(W2); title("from Euclidean distance")
#' par(opar)
#' }
#'
#' @author Kisung You
#' @references Floyd, R.W. (1962) \emph{Algorithm 97: Shortest Path}. Commincations of the ACMS, Vol.5(6):345.
#' @rdname aux_shortestpath
#' @export
aux.shortestpath <- function(dist){
  # class determination
  distnaive = as.matrix(dist)
  if ((nrow(distnaive)!=ncol(distnaive))||(!isSymmetric(distnaive))){
    stop("* aux.shortestpath : input 'dist' should be either (n*n) matrix or 'dist' class object.")
  }
  # consider logical input
  if (any(is.logical(distnaive))){
    distnaive = distnaive*1
  }
  # set as -Inf for 0 values
  mepsil  = .Machine$double.eps
  distnaive[which(distnaive<5*mepsil)] = -Inf
  distgeo   = aux_shortestpath(distnaive)
  return(distgeo)
}
