#' Manifold-Adaptive Dimension Estimation
#'
#' \code{do.made} first aims at finding local dimesion estimates using nearest neighbor techniques based on
#' the first-order approximation of the probability mass function and then combines them to get a single global estimate. Due to the rate of convergence of such
#' estimate to be independent of assumed dimensionality, authors claim this method to be
#' \emph{manifold-adaptive}.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param k size of neighborhood for analysis.
#' @param maxdim maximum possible dimension allowed for the algorithm to investigate.
#' @param combine method to aggregate local estimates for a single global estimate.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{estimated global intrinsic dimension.}
#' \item{locdim}{a length-\eqn{n} vector estimated dimension at each point.}
#' }
#'
#' @examples
#' \dontrun{
#' ## create a data set of intrinsic dimension 2.
#' X = aux.gensamples(dname="swiss")
#'
#' ## compare effect of 3 combining scheme
#' out1 = est.made(X, combine="mean")
#' out2 = est.made(X, combine="median")
#' out3 = est.made(X, combine="vote")
#'
#' ## print the results
#' sprintf("* est.made : estimated dimension with 'mean'   method is %d.",out1$estdim)
#' sprintf("* est.made : estimated dimension with 'median' method is %d.",out2$estdim)
#' sprintf("* est.made : estimated dimension with 'vote'   method is %d.",out3$estdim)
#' }
#'
#' @references
#' \insertRef{farahmand_manifold-adaptive_2007}{Rdimtools}
#'
#' @author Kisung You
#' @export
est.made <- function(X, k=round(sqrt(ncol(X))), maxdim=min(ncol(X),15), combine=c("mean","median","vote")){
  ##########################################################################
  ## Preprocessing and Default Parameter
  aux.typecheck(X)
  n = nrow(X)
  D = as.matrix(dist(X))
  diag(D) = Inf

  ##########################################################################
  ## Computation
  #   1. sort D
  Dsort = apply(D, 2, sort)
  #   2. Rk and Rk2
  Rk    = as.vector(Dsort[k,])
  Rk2   = as.vector(Dsort[ceiling(k/2),])
  #   3. compute a pointwise estimate
  dhat  = log(2)/(log(Rk/Rk2))
  #   4. combine
  combine = match.arg(combine)
  estdim  = switch(combine,
                   mean   = max(round(mean(pmin(dhat,maxdim))),1),
                   median = max(round(median(pmin(dhat,maxdim))),1),
                   vote   = max(as.integer(names(which.max(table(round(dhat)))[1])),1)
                   )

  ##########################################################################
  ## Report
  result = list()
  result$estdim = estdim
  result$locdim = dhat
  return(result)
}
