#' Maximum Likelihood Esimation with Poisson Process
#'
#' Assuming the density in a hypersphere is constant, authors proposed to build
#' a likelihood structure based on modeling local spread of information via Poisson Process.
#' \code{est.mle1} requires two parameters that model the reasonable range of neighborhood size
#' to reflect inhomogeneity of distribution across data points.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param k1 minimum neighborhood size, larger than 1.
#' @param k2 maximum neighborhood size, smaller than \eqn{n}.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{estimated intrinsic dimension.}
#' }
#'
#' @examples
#' \donttest{
#' ## create example data sets with intrinsic dimension 2
#' X1 = aux.gensamples(dname="swiss")
#' X2 = aux.gensamples(dname="ribbon")
#' X3 = aux.gensamples(dname="saddle")
#'
#' ## acquire an estimate for intrinsic dimension
#' out1 = est.mle1(X1)
#' out2 = est.mle1(X2)
#' out3 = est.mle1(X3)
#'
#' ## print the estimates
#' line1 = paste0("* est.mle1 : 'swiss'  estiamte is ",round(out1$estdim,2))
#' line2 = paste0("* est.mle1 : 'ribbon' estiamte is ",round(out2$estdim,2))
#' line3 = paste0("* est.mle1 : 'saddle' estiamte is ",round(out3$estdim,2))
#' cat(paste0(line1,"\n",line2,"\n",line3))
#' }
#'
#' @references
#' \insertRef{levina_maximum_2005}{Rdimtools}
#'
#' @rdname estimate_mle1
#' @author Kisung You
#' @export
est.mle1 <- function(X, k1=10, k2=20){
  ##########################################################################
  ## Preprocessing and Default Parameter
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  d = stats::dist(X)
  D = as.matrix(d)
  for (i in 1:n){
    D[i,i]=Inf
  }
  k1 = as.integer(k1)
  k2 = as.integer(k2)
  if ((length(k1)>=2)||(length(k2)>=2)||(k2>=n)||(k1<=1)||(k1>=k2)){
    stop("* est.mle1 : 'k1' and 'k2' should be suitably selected.")
  }

  ##########################################################################
  ## Comutation
  #   1. get knn distance information
  matK  = est.nearkdist(D, k2)
  #   2. iterate over k
  mkhat = rep(0,(k2-k1+1))
  for (k in k1:k2){
    mkhat[k-k1+1] = mean(1/((log(matK[,k]) - rowMeans(log(matK)[, 1:(k-1)]))))
  }
  #   3. take average
  estdim = sum(mkhat)/(k2-k1+1)

  ##########################################################################
  ## Return Result
  result = list()
  result$estdim = estdim
  return(result)
}


#   -----------------------------------------------------------------------
#' @keywords internal
#' @noRd
est.nearkdist <- function(D, k){
  n = nrow(D)
  output = array(0,c(n,k))
  for (i in 1:n){
    output[i,] = sort(as.vector(D[i,]))[1:k]
  }
  return(output)
}
