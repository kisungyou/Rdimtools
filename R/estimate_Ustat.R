#' ID Estimation with Convergence Rate of U-statistic on Manifold
#'
#' \eqn{U}-statistic is built upon theoretical arguments with the language of
#' smooth manifold. The convergence rate of the statistic is achieved as a proxy
#' for the estimated dimension by, at least partially, considering
#' the scale and influence of extrinsic curvature. The method returns \emph{integer} valued
#' estimate in that there is no need for rounding the result for practical usage.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param maxdim maximum possible dimension allowed for the algorithm to investigate.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{estimated intrinsic dimension.}
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
#' out1 = est.Ustat(X1)
#' out2 = est.Ustat(X2)
#' out3 = est.Ustat(X3)
#'
#' ## print the results
#' line1 = paste0("* est.Ustat : 'swiss'  gives ",round(out1$estdim,2))
#' line2 = paste0("* est.Ustat : 'ribbon' gives ",round(out2$estdim,2))
#' line3 = paste0("* est.Ustat : 'saddle' gives ",round(out3$estdim,2))
#' cat(paste0(line1,"\n",line2,"\n",line3))
#' }
#'
#' @references
#' \insertRef{hein_intrinsic_2005}{Rdimtools}
#'
#' @rdname estimate_Ustat
#' @author Kisung You
#' @export
est.Ustat <- function(X, maxdim=min(ncol(X),15)){
  ##########################################################################
  ## Preprocessing and Default Parameter
  aux.typecheck(X)
  N = nrow(X)
  D = as.matrix(dist(X))

  ##########################################################################
  ## Preliminary Computations
  #   1. hlN : d(X_i)
  hlN = sum(apply(D,1,sort)[2,])/N
  #   2. kernel function (cdim : current dimensionality)
  kh <- function(x, h, cdim) {
    ifelse((1 - (x/h)^2) > 0, 1 - (x/h)^2, 0)/(h^cdim)
  }

  ##########################################################################
  ## Main Computation : iterate over dimension
  slopes = rep(0,maxdim) # following the notation for dimension as 'l'
  for (l in 1:maxdim){
    Ustats = rep(0,5)

    # 1. manually for the r=1
    Ustats[1] = mean(kh(as.numeric(as.dist(D)), 1, l))
    # 2. for r=2 to 5, use 'aux.randpartition' and 'combn'
    for (r in 2:5){
      Ustats[r] = est.Ustat.singler(D, l, r, kh, hlN)
    }
    # 3. regression
    ns = round(N/(1:5))

    x  = log(hlN * ((N/ns) * (log(ns)/log(N)))^l)
    y  = log(Ustats)
    finiteid = intersect(which(is.finite(x)), which(is.finite(y)))

    x = x[finiteid]
    y = y[finiteid]

    slopes[l] = abs(sum(lm(y~x+1,weights=(1/(1:5))[finiteid])$coef[2]))
  }

  ##########################################################################
  ## Wrap Up
  result = list()
  result$estdim = max(which.min(slopes), 1)
  return(result)
}



## single iterate with l, r
#' @keywords internal
#' @noRd
est.Ustat.singler <- function(D, l, r, func, hlN){
  N       = nrow(D)
  indices = aux.randpartition(N,r) # length-r list of indices
  n       = round(N/r)
  hln     = hlN * ((N/n) * (log(n)/log(N)))^l

  allcase = combn(1:r,2)
  ncases  = ncol(allcase)
  Ustats  = rep(0,ncases)

  for (iter in 1:ncases){
    id1 = as.integer(allcase[1,iter]) # index of randomly partitioned indices
    id2 = as.integer(allcase[2,iter])

    idx1 = indices[[id1]]
    idx2 = indices[[id2]]

    tmpD = D[idx1,idx2]
    if (id1==id2){
      Ustats[iter] = mean(func(as.numeric(as.dist(tmpD)), hln, l))
    } else {
      Ustats[iter] = mean(func(as.numeric(tmpD), hln, l))
    }
  }
  return(mean(Ustats))
}
