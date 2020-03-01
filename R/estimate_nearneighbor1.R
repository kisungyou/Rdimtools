#' Intrinsic Dimension Estimation with Near-Neighbor Information
#'
#' Based on an assumption of data points being locally uniformly distributed,
#' \code{est.nearneighbor1} estimates the intrinsic dimension based on the
#' local distance information in an iterative manner.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param K maximum neighborhood size, smaller than \eqn{p}.
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
#' output = est.nearneighbor1(X)
#' sprintf("* est.nearneighbor1 : estimated dimension is %.2f.",output$estdim)
#' }
#'
#' @references
#' \insertRef{pettis_intrinsic_1979}{Rdimtools}
#'
#' @author Kisung You
#' @export
est.nearneighbor1 <- function(X, K=max(2,round(ncol(X)/5))){
  ##########################################################################
  ## Preprocessing and Default Parameter
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  d = stats::dist(X)
  D = as.matrix(d)
  stopthr = 0.01
  if ((K>=ncol(X))||(K<1)){
    stop("* est.nearneighbor1 : 'K' should be in (1,ncol(X)).")
  }
  maxiter = 100

  ##########################################################################
  ## Preliminary Computation
  nbdmat = matrix(0, ncol = K, nrow = n)
  trcval = nbdmat
  for (i in 1:n) {
    nbdmat[i, ] = order(D[i, ])[2:(K + 1)]
    trcval[i, ] = D[i, nbdmat[i, ]]
  }
  trcval = trcval[, 1:K]
  barRk  = colMeans(trcval)

  ##########################################################################
  ## Main Computation
  #   1. initialization
  vecK    = seq(from=1,to=K,length.out=K)
  dest    = 1/sum(stats::coefficients(stats::lm(log(barRk)~log(vecK)))[2])
  destvec = dest

  #   2. main iteration
  for (k in 1:maxiter) {
    logGkd = est.nni.logGkd(K, dest) # update
    dest   = 1/sum(stats::coefficients(stats::lm((log(barRk)+logGkd)~log(vecK)))[2])
    if (abs(destvec[length(destvec)]-dest) < stopthr) {
      break
    }
    destvec <- c(destvec, dest)      # update
  }

  #####################################
  # STEP 4 : return the result
  output = list()
  output$estdim = dest
  return(output)
}

#' @keywords internal
#' @noRd
est.nni.logGkd <- function(K, d){
  vecK  = 1:K
  term1 = (d-1)/(2*(d^2)*vecK)
  term2 = ((d-1)*(d-2))/(12*(d^3)*(vecK^2))
  term3 = ((d-1)^2)/(12*(d^4)*(vecK^3))
  term4 = ((d-1)*(d-2)*((d^2)+(3*d)-3))/(120*(d^5)*(vecK^4))
  output = term1+term2-term3-term4
  return(output)
}
