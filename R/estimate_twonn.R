#' Intrinsic Dimension Estimation by a Minimal Neighborhood Information
#'
#' Unlike many intrinsic dimension (ID) estimation methods, \code{est.twonn} only requires
#' two nearest datapoints from a target point and their distances. This extremely minimal approach
#' is claimed to redue the effects of curvature and density variation across different locations
#' in an underlying manifold.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
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
#' out1 = est.twonn(X1)
#' out2 = est.twonn(X2)
#' out3 = est.twonn(X3)
#'
#' ## print the results
#' sprintf("* est.twonn : estimated dimension for 'swiss'  data is %.2f.",out1$estdim)
#' sprintf("* est.twonn : estimated dimension for 'ribbon' data is %.2f.",out2$estdim)
#' sprintf("* est.twonn : estimated dimension for 'saddle' data is %.2f.",out3$estdim)
#' }
#'
#' @references
#' \insertRef{facco_estimating_2017}{Rdimtools}
#'
#' @rdname estimate_twonn
#' @author Kisung You
#' @export
est.twonn <- function(X){
  ##########################################################################
  ## Preprocessing and Default Parameter
  aux.typecheck(X)
  n = nrow(X)
  D = as.matrix(dist(X))
  diag(D) = Inf

  ##########################################################################
  ## Computation
  #   1. compute the ratio
  mu = rep(0,n)
  for (i in 1:n){
    tgt = sort(as.vector(D[i,]))[1:2]
    mu[i] = tgt[2]/tgt[1]
  }
  #   2. transform
  mu.sorted = sort(mu)
  empiF     = (1:(n-1))/n
  #   3. let's do the linear regression
  x = log(mu.sorted)[1:(n-1)]
  y = -log(1.0-empiF)

  ##########################################################################
  ## Return the results
  result = list()
  result$estdim = sum(coefficients(lm(y~x-1))[1])
  return(result)
}
