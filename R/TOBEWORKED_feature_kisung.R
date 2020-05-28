#' Kisung's
#'
#' idea
#' 1) compute distance correlation
#' 2) quantile and threshold (Aij = 1 if DCorr(i,j) >= 'x' quantile)
#' 3) shortest path
#' 4) PAM -> partition id -> use as feature indices
#'
#'
#' @keywords internal
do.kisung <- function(X, ndim=2, qthr=0.75){
  # library(T4CPR)
  # X   = as.matrix(iris[,1:4])
  # drX = T4CPR::estcorr(X, method="Distance")$R
  # vecdX = as.vector(drX[upper.tri(drX)])
  #
  # thr = as.double(quantile(vecdX, probs=0.25))
  # A   = array(0,c(nrow(drX),ncol(drX)))
  # A[(drX>thr)] = 1; diag(A) = 0
  # spA = aux.shortestpath(A)
  # ddA = as.dist(spA)
  #
  # cluster::pam(ddA, k=2)$id.med

  # library(PNAS)
  # library(T4CPR)
  # X   = as.matrix(iris[,1:4])
  # drX = T4CPR::estcorr(X, method="Distance")$R; diag(drX)=0
  #
  # thr = 0.85
  # A = array(0,c(4,4))
  # A[(drX > thr)] = 1
  # L = diag(rowSums(A))-A
  # eigen(L)
  # plot(eigen(L)$vectors[,1:2])
  # plot(eigen(L)$vectors[,3:4])
  #
  # g0 = graph_from_adjacency_matrix(array(0,c(4,4)),mode = "undirected")
}
