#' Intrinsic Dimension Estimation using Packing Numbers
#'
#' Instead of covering numbers which are expensive to compute in many fractal-based methods,
#' \code{est.packing} exploits packing numbers as a proxy to describe spatial density. Since
#' it involves random permutation of the dataset at each iteration, every run might have
#' different results.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param eps small positive number for stopping threshold.
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
#' ## try different eps values
#' out1 = est.packing(X, eps=0.1)
#' out2 = est.packing(X, eps=0.01)
#' out3 = est.packing(X, eps=0.001)
#'
#' ## print the results
#' line1 = paste0("* est.packing : eps=0.1   gives ",round(out1$estdim,2))
#' line2 = paste0("* est.packing : eps=0.01  gives ",round(out2$estdim,2))
#' line3 = paste0("* est.packing : eps=0.001 gives ",round(out3$estdim,2))
#' cat(paste0(line1,"\n",line2,"\n",line3))
#' }
#'
#' @references
#' \insertRef{kegl_intrinsic_2002}{Rdimtools}
#'
#' @rdname estimate_packing
#' @author Kisung You
#' @export
est.packing <- function(X, eps=0.01){
  ##########################################################################
  ## Preprocessing and Default Parameter
  aux.typecheck(X)
  n  = nrow(X)
  D  = as.matrix(dist(X))
  qD = as.numeric(stats::quantile(D[upper.tri(D)], seq(from=0,to=1,length.out=11))) # quantile and triangular part

  rvec = qD[2:3] # 10% and 20% quantiles

  ##########################################################################
  ## Let's follow the pseudocode in the paper
  flag  = FALSE
  l     = 0 # I will iterate over
  matL  = c()
  while (flag==FALSE){
    l  = l+1                             # update l for 1 to \infty
    X  = X[sample(1:n,n,replace=FALSE),] # permute Sn randomly
    Lk = c()
    for (k in 1:2){
      rk = rvec[k]
      C  = sample(1:n, 1)               # we can't run as pseudocode directly.
      for (i in 2:n){
        idC = c()
        for (j in C){
          idC = c(idC, which(D[j,]<rk))
        }
        idC = unique(idC)
        if (length(idC)<n){
          C = rbind(C, setdiff(1:n, idC)[1])
        } else {
          break
        }
      }
      Lk = c(Lk, log(length(C)))
    }
    matL  = cbind(matL, Lk)
    Dpack = -(Lk[2]-Lk[1])/(log(rvec[2])-log(rvec[1]))

    cond1 = (l>10)
    cond2 = (1.65*sqrt(sum(apply(matL,1,stats::var)))/sqrt(l*(log(rvec[2])-log(rvec[1]))) < Dpack*(1-eps)/2)
    flag  = (cond1&&cond2)
  }

  result = list()
  result$estdim = max(Dpack,1)
  return(result)
}
