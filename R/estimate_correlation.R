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
#' Even though we could use arbitrary \code{cut} to compute estimated dimension, it is also possible to
#' use visual inspection. According to the theory, if the function returns an \code{output}, we can plot
#' \code{plot(log(output$r), log(output$Cr))} and use the linear slope in the middle as desired dimension of data.
#'
#' @section Automatic choice of \eqn{r}:
#' The least value for radius \eqn{r} must have non-degenerate counts, while the maximal value should be the
#' maximum distance among all pairs of data points across all coordinates. \code{nlevel} controls the number of interim points
#' in a log-equidistant manner.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param nlevel the number of \code{r} (radius) to be tested.
#' @param cut a vector of ratios for computing estimated dimension in \eqn{(0,1)}.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{estimated dimension using \code{cut} values.}
#' \item{r}{a vector of radius used.}
#' \item{Cr}{a vector of \eqn{C(r)} as decribed above.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate three different dataset
#' X1 = aux.gensamples(dname="swiss")
#' X2 = aux.gensamples(dname="ribbon")
#' X3 = aux.gensamples(dname="twinpeaks")
#'
#' ## compute
#' out1 = est.correlation(X1)
#' out2 = est.correlation(X2)
#' out3 = est.correlation(X3)
#'
#' ## visually verify : all should have approximate slope of 2.
#' par(mfrow=c(1,3))
#' plot(log(out1$r), log(out1$Cr), main="swiss roll")
#' plot(log(out2$r), log(out2$Cr), main="ribbon")
#' plot(log(out3$r), log(out3$Cr), main="twinpeaks")
#' }
#'
#'
#' @references
#' \insertRef{grassberger_measuring_1983}{Rdimtools}
#'
#' @author Kisung You
#' @seealso \code{\link{est.boxcount}}
#' @rdname estimate_correlation
#' @export
est.correlation <- function(X,nlevel=50,cut=c(0.1,0.9)){
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
  conversion = approx_correlation(vecr, vecCm, nlevel)
  vecr   = conversion$r
  vecCm  = conversion$Cr

  #   5-2. estimated dimension
  if ((!is.vector(cut))||(length(cut)!=2)||(any(cut<=0))||(any(cut>=1))){
    stop("* est.correlation : 'cut' should be a vector of length 2.")
  }
  idxtrim = fractal_trimmer(vecCm, cut)
  vecrt   = vecr[idxtrim]  # trimmed vecr
  vecCt   = vecCm[idxtrim] # trimmed vecCm
  nnn     = length(vecrt)  # number of saved ones



  estdim   = sum(coef(lm(log(vecCt)~log(vecrt)))[2])


  #   5-3. show and return
  output = list()
  output$estdim = estdim
  output$r  = vecr
  output$Cr = vecCm
  return(output)
}


#' @keywords internal
#' @noRd
approx_correlation <- function(r, Cr, nlevel){
  x = log(r)
  y = log(Cr)

  interp = stats::approx(x,y,n=nlevel)
  xnew   = interp$x
  ynew   = interp$y

  output = list()
  output$r  = exp(xnew)
  output$Cr = exp(ynew)
  return(output)
}
