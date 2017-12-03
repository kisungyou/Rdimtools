#' Correlation Dimension
#'
#' Correlation dimension is a measure of determining the dimension of a given set. It is
#' often referred to as a type of fractal dimension.  Its mechanism is somewhat similar to
#' that of box-counting dimension, but has the advantage of being intuitive as well as
#' efficient in terms of computation with some robustness contingent on the lack of availability for large dataset.
#' \deqn{dim(S) = \lim \frac{\log C(r)}{\log r}} as \eqn{r\rightarrow 0}, where
#' \eqn{C(r)=\lim (2/(N-1)*N)\sum_i^N \sum_{j=i+1}^N I(\|x_i-x_j\|\le r)}.
#'
#' @section Determining the dimension:
#' In this version, no automated method is included due to lots of possibilities for adopting
#' any types of change point detection methods. Instead, we recommend using \emph{visual}
#' identification by looking at the slope of the linear part in the middle using either
#' \code{show=TRUE} flag or plotting as \code{plot(log(output$r),log(output$Cr))}.
#'
#' @param X an \code{(n-by-p)} matrix or data frame whose rows are observations.
#' @param nlevel the number of \code{r} (radius) to be tested.
#' @param show a logical; \code{FALSE} for not showing visually, \code{TRUE} otherwise.
#'
#' @return a named data frame containing \describe{
#' \item{r}{a vector of radius used.}
#' \item{Cr}{a vector of \eqn{C(r)} as decribed above.}
#' }
#' @examples
#' ## generate three different dataset
#' X1 = aux.gensamples(dname="swiss")
#' X2 = aux.gensamples(dname="ribbon")
#' X3 = aux.gensamples(dname="twinpeaks")
#'
#' ## visually verify : all should have approximate slope of 2.
#' par(mfrow=c(1,3))
#' est.correlation(X1,show=TRUE)
#' est.correlation(X2,show=TRUE)
#' est.correlation(X3,show=TRUE)
#'
#'
#' @references Grassberger, P. and Procaccia, I. (1983) \emph{Measuring the strangeness of strange attractors}. Physica D9:189-208.
#' @author Kisung You
#' @seealso \code{\link{est.boxcount}}
#' @rdname estimate_correlation
#' @export

est.correlation <- function(X,nlevel=50,show=FALSE){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  X = as.matrix(X)
  n = nrow(X)
  d = ncol(X)
  if ((!is.numeric(nlevel))||is.infinite(nlevel)||is.na(nlevel)||(nlevel<2)){
    stop("* est.boxcount : 'nlevel' should be a positive integer bigger than 1.")
  }
  nlevel = as.integer(nlevel)

  # 2. min/max for each dimension
  Dmat   = (dist(X))
  rstart = max(Dmat)+0.05
  rmin   = max(min(Dmat),0.01)

  # 3. setting for possible computations
  vecr    = exp(seq(from=log(rstart),to=log(rmin),length.out=nlevel))
  rlength = length(vecr)
  vecCm   = as.vector(array(0,c(1,rlength)))

  # 4. main iteration
  n2 = 2/(n*(n-1))
  for (i in 1:rlength){
    # 4-1. target
    tgtr = vecr[i]
    # 4-2. count
    vecCm[i] = length(which(Dmat<=tgtr))*n2
  }

  # 5. return output
  #   5-1. ratio
  ratio = log(vecCm)/log(vecr)
  output = data.frame(r=vecr,Cr=vecCm)

  #   5-3. show and return
  result = output
  if (show){
    plot(log(output$r),log(output$Cr),
         main="correlation::linear slope",
         xlab = "log(r)", ylab = "log(C(r))",type="b")
  }
  return(result)
}
