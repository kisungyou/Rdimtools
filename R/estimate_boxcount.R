#' Box-counting dimension
#'
#' Box-counting dimension, also known as Minkowski-Bouligand dimension, is a popular way of figuring out
#' the fractal dimension of a set in a Euclidean space. Its idea is to measure the number of boxes
#' required to cover the set repeatedly by decreasing the length of each side of a box. It is defined as
#' \deqn{dim(S) = \lim \frac{\log N(r)}{\log (1/r)}} as \eqn{r\rightarrow 0}, where \eqn{N(r)} is
#' the number of boxes counted to cover a given set for each corresponding \eqn{r}.
#'
#' @section Determining the dimension:
#' In this version, no automated method is included due to lots of possibilities for adopting
#' any types of change point detection methods. Instead, we recommend using \emph{visual}
#' identification by looking at the slope of the linear part in the middle using either
#' \code{show=TRUE} flag or plotting as \code{plot(log(1/output$r),log(output$Nr))}.
#'
#' @param X an \code{(n-by-p)} matrix or data frame whose rows are observations.
#' @param nlevel the number of \code{r} (radius) to be tested.
#' @param show a logical; \code{FALSE} for not showing visually, \code{TRUE} otherwise.
#'
#' @return a named data frame containing \describe{
#' \item{r}{a vector of radius used.}
#' \item{Nr}{a vector of boxes counted for each corresponding \code{r}.}
#' }
#' @examples
#' ## generate three different dataset
#' X1 = aux.gensamples(dname="swiss")
#' X2 = aux.gensamples(dname="ribbon")
#' X3 = aux.gensamples(dname="twinpeaks")
#'
#' ## visually verify : all should have approximate slope of 2.
#' par(mfrow=c(1,3))
#' est.boxcount(X1,show=TRUE)
#' est.boxcount(X2,show=TRUE)
#' est.boxcount(X3,show=TRUE)
#'
#' @references Ott, E. (1988) \emph{Chaos in Dynamical Systems}. Cambridge University Press.
#' @author Kisung You
#' @seealso \code{\link{est.correlation}}
#' @rdname estimate_boxcount
#' @export
est.boxcount <- function(X,nlevel=50,show=FALSE){
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
  Is = aux_minmax(X,0.1)
  rstart = ((min(Is[2,]-Is[1,])) + (max(Is[2,]-Is[1,])))/2

  # 3. setting for possible computations
  maxlength = floor(1-log2(max(1e-15,1e-10/rstart)))
  maxlength = min(nlevel,maxlength)

  # 4. main iteration
  vecr = as.vector(array(0,c(1,maxlength)))
  vecN = as.vector(array(0,c(1,maxlength)))
  Imin = as.vector(Is[1,])
  tX   = t(X)
  for (i in 1:maxlength){
    currentr = rstart*(0.75^(i-1))
    counted  = round(methods_boxcount(tX,Imin,currentr))
    vecr[i] = currentr
    vecN[i] = sum(!duplicated(counted))
  }

  # 5. return output
  #   5-0. trimming of Nr max
  if (length(which(vecN==max(vecN)))>2){
    minidx = min(which(vecN==max(vecN)))
    vecr = vecr[1:(minidx+1)]
    vecN = vecN[1:(minidx+1)]
  }
  #   5-1. ratio
  ratio = log(vecN)/log(1/vecr)
  output = data.frame(r=vecr, Nr=vecN)

  #   5-3. show and return
  result = output
  if (show){
    plot(log(1/output$r),log(output$Nr),
         main="boxcount:: use linear slope",
         xlab = "log(1/r)", ylab = "log(N(r))",type="b")
  }
  return(result)
}
