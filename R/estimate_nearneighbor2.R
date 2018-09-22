#' Near-Neighbor Information with Bias Correction
#'
#' Though similar to \code{\link{est.nearneighbor1}}, authors of the reference
#' argued that there exists innate bias in the method and proposed a non-iterative algorithm
#' to reflect local distance information under a range of neighborhood sizes.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param kmin minimum neighborhood size, larger than 1.
#' @param kmax maximum neighborhood size, smaller than \eqn{p}.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{estimated intrinsic dimension.}
#' }
#'
#' @examples
#' \dontrun{
#' ## create an example data with intrinsic dimension 2
#' X = cbind(aux.gensamples(dname="swiss"),aux.gensamples(dname="swiss"))
#'
#' ## acquire an estimate for intrinsic dimension
#' output = est.nearneighbor2(X)
#' sprintf("* est.nearneighbor2 : estimated dimension is %.2f.",output$estdim)
#' }
#'
#' @references
#' \insertRef{verveer_evaluation_1995}{Rdimtools}
#'
#' @export
est.nearneighbor2 <- function(X, kmin=2, kmax=max(3, round(ncol(X)/2))){
  ##########################################################################
  ## Preprocessing and Default Parameter
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  d = stats::dist(X)
  D = as.matrix(d)
  stopthr = 0.01
  if ((kmin<2)||(kmin>=kmax)||(kmax>=ncol(X))){
    stop("* est.nearneighbor2 : selection of 'kmin' and 'kmax' are incorrect.")
  }
  maxiter = 100

  ##########################################################################
  ## Preliminary Computation
  nbdmat = matrix(0, ncol = kmax, nrow = n)
  trcval = nbdmat
  for (i in 1:n) {
    nbdmat[i, ] = order(D[i, ])[2:(kmax + 1)]
    trcval[i, ] = D[i, nbdmat[i, ]]
  }
  trcval = trcval[, 1:kmax]
  barRk  = colMeans(trcval)

  ##########################################################################
  ## main computation
  #   term1 : to be inverted
  term1 = 0
  for (k in kmin:(kmax-1)){
    term1 = ((barRk[k+1]-barRk[k])^2) + term1
  }
  #   term2 : scaled difference
  term2 = 0
  for (k in kmin:(kmax-1)){
    term2 = term2 + (((barRk[k+1]-barRk[k])*barRk[k])/k)
  }

  #####################################
  # return the result
  output = list()
  output$estdim = (term2/term1)
  return(output)
}
