#' MINDml
#'
#' It is a minimum neighbor distance estimator of the intrinsic dimension based on Maximum Likelihood principle.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param k the neighborhood size for defining locality.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{the global estimated dimension.}
#' }
#'
#' @examples
#' ## create 3 datasets of intrinsic dimension 2.
#' set.seed(100)
#' X1 = aux.gensamples(dname="swiss")
#' X2 = aux.gensamples(dname="ribbon")
#' X3 = aux.gensamples(dname="saddle")
#'
#' ## acquire an estimate for intrinsic dimension
#' out1 = est.mindml(X1, k=10)
#' out2 = est.mindml(X2, k=10)
#' out3 = est.mindml(X3, k=10)
#'
#' ## print the results
#' line1 = paste0("* est.mindml : 'swiss'  estiamte is ",round(out1$estdim,2))
#' line2 = paste0("* est.mindml : 'ribbon' estiamte is ",round(out2$estdim,2))
#' line3 = paste0("* est.mindml : 'saddle' estiamte is ",round(out3$estdim,2))
#' cat(paste0(line1,"\n",line2,"\n",line3))
#'
#' @references
#' \insertRef{lombardi_minimum_2011}{Rdimtools}
#'
#' @seealso \code{\link{est.mindkl}}
#'
#' @rdname estimate_mindml
#' @author Kisung You
#' @export
est.mindml <- function(X, k=5){
  ##########################################################################
  ## preprocessing
  aux.typecheck(X)
  N = nrow(X)
  D = ncol(X)
  k = round(k)

  ##########################################################################
  ## computation taken from DANCo
  dX = as.matrix(stats::dist(X))
  vec.topK1 <- list()
  for (n in 1:N){
    tgt = as.vector(dX[n,])
    vec.topK1[[n]] = order(tgt)[2:(k+2)]
  }
  data.dML = danco_part1(dX, vec.topK1, N, D, k)

  ##########################################################################
  ## Return the results
  result = list()
  result$estdim = data.dML
  return(result)
}
