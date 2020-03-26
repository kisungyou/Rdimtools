#' MiNDkl
#'
#' It is a minimum neighbor distance estimator of the intrinsic dimension based on Kullback Leibler divergence estimator.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param k the neighborhood size for defining locality.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{the global estimated dimension.}
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
#' out1 = est.mindkl(X1, k=10)
#' out2 = est.mindkl(X2, k=10)
#' out3 = est.mindkl(X3, k=10)
#'
#' ## print the results
#' sprintf("* est.mindkl : estimated dimension for 'swiss'  data is %.2f.",out1$estdim)
#' sprintf("* est.mindkl : estimated dimension for 'ribbon' data is %.2f.",out2$estdim)
#' sprintf("* est.mindkl : estimated dimension for 'saddle' data is %.2f.",out3$estdim)
#' }
#'
#' @references
#' \insertRef{lombardi_minimum_2011}{Rdimtools}
#'
#' @seealso \code{\link{est.mindml}}
#'
#' @rdname estimate_mindkl
#' @author Kisung You
#' @export
est.mindkl <- function(X, k=5){
  ##########################################################################
  ## preprocessing
  aux.typecheck(X)
  N = nrow(X)
  D = ncol(X)
  k = round(k)

  ##########################################################################
  ## preliminary computation
  dX = as.matrix(stats::dist(X))
  vec.topK1 <- list()
  for (n in 1:N){
    tgt = as.vector(dX[n,])
    vec.topK1[[n]] = order(tgt)[2:(k+2)]
  }
  rho.data = mindkl_distance(dX, vec.topK1, N, D, k)

  ##########################################################################
  ## iteration
  cost.d <- rep(0,D)
  for (d in 1:D){
    # generate data
    Y  = matrix(stats::rnorm(N*d),ncol=d)
    Y  = (Y/base::sqrt(base::rowSums(Y^2)))*(stats::runif(1, min=1e-8, max=1-(1e-8))^(1/d))

    # extract neighborhood information & rho vector
    dY = as.matrix(stats::dist(Y))
    vec.topK1 <- list()
    for (n in 1:N){
      tgt = as.vector(dY[n,])
      vec.topK1[[n]] = order(tgt)[2:(k+2)]
    }
    rho.yd = mindkl_distance(dY, vec.topK1, N, D, k)

    # record the cost
    cost.d[d] = (log(N)-log(N-1)) + (sum(log(rho.data)-log(rho.yd))/N)
  }

  ##########################################################################
  ## report
  idmin = which.min(cost.d)
  if (length(idmin) > 1){
    d.fin = sample(idmin, 1)
  } else {
    d.fin = idmin
  }
  result = list()
  result$estdim = d.fin
  return(result)
}


# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
mindkl_distance <- function(dmat, vec.topK1, N, D, k){
  vec.rho = rep(0,N)
  for (n in 1:N){
    tgt = dmat[n,vec.topK1[[n]]]
    vec.rho[n] = min(tgt/(max(tgt)))
  }
  heyo = as.matrix(stats::dist(matrix(vec.rho, ncol=1)))
  output = rep(0,N)
  for (n in 1:N){
    tgt = base::sort(as.vector(heyo[n,]))
    output[n] = tgt[2]
  }
  return(output)
}
