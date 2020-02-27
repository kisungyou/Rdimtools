#' MINDml
#'
#'
#' @examples
#' \dontrun{
#' ## create 3 datasets of intrinsic dimension 2.
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
#' sprintf("* est.mindml : estimated dimension for 'swiss'  data is %.2f.",out1$estdim)
#' sprintf("* est.mindml : estimated dimension for 'ribbon' data is %.2f.",out2$estdim)
#' sprintf("* est.mindml : estimated dimension for 'saddle' data is %.2f.",out3$estdim)
#' }
#'
#' @references
#' \insertRef{lombardi_minimum_2011}{Rdimtools}
#'
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
