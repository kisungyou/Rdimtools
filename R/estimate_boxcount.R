#' Box-counting Dimension
#'
#' Box-counting dimension, also known as Minkowski-Bouligand dimension, is a popular way of figuring out
#' the fractal dimension of a set in a Euclidean space. Its idea is to measure the number of boxes
#' required to cover the set repeatedly by decreasing the length of each side of a box. It is defined as
#' \deqn{dim(S) = \lim \frac{\log N(r)}{\log (1/r)}} as \eqn{r\rightarrow 0}, where \eqn{N(r)} is
#' the number of boxes counted to cover a given set for each corresponding \eqn{r}.
#'
#' @section Determining the dimension:
#' Even though we could use arbitrary \code{cut} to compute estimated dimension, it is also possible to
#' use visual inspection. According to the theory, if the function returns an \code{output}, we can plot
#' \code{plot(log(1/output$r),log(output$Nr))} and use the linear slope in the middle as desired dimension of data.
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
#' \item{estdim}{estimated dimension using \code{cut} ratios.}
#' \item{r}{a vector of radius used.}
#' \item{Nr}{a vector of boxes counted for each corresponding \code{r}.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate three different dataset
#' X1 = aux.gensamples(dname="swiss")
#' X2 = aux.gensamples(dname="ribbon")
#' X3 = aux.gensamples(dname="twinpeaks")
#'
#' ## compute boxcount dimension
#' out1 = est.boxcount(X1)
#' out2 = est.boxcount(X2)
#' out3 = est.boxcount(X3)
#'
#' ## visually verify : all should have approximate slope of 2.
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(log(1/out1$r), log(out1$Nr), main="swiss roll")
#' plot(log(1/out2$r), log(out2$Nr), main="ribbon")
#' plot(log(1/out3$r), log(out3$Nr), main="twinpeaks")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{hentschel_infinite_1983}{Rdimtools}
#'
#' \insertRef{ott_chaos_2002}{Rdimtools}
#'
#' @author Kisung You
#' @seealso \code{\link{est.correlation}}
#' @rdname estimate_boxcount
#' @export
est.boxcount <- function(X,nlevel=50,cut=c(0.1,0.9)){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  X = as.matrix(X)
  n = nrow(X)
  d = ncol(X)
  nlevel = as.integer(nlevel)
  if ((!is.numeric(nlevel))||is.infinite(nlevel)||is.na(nlevel)||(nlevel<2)){
    stop("* est.boxcount : 'nlevel' should be a positive integer bigger than 1.")
  }

  # 2. min/max for each dimension
  Is     = aux_minmax(X,0.1)
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
  #   5-1. conversion
  conversion = approx_boxcount(vecr, vecN, nlevel)
  vecr   = conversion$r
  vecN   = conversion$Nr

  #   5-2. estimated dimension with trimming
  if ((!is.vector(cut))||(length(cut)!=2)||(any(cut<=0))||(any(cut>=1))){
    stop("* est.boxcount : 'cut' should be a vector of length 2.")
  }

  idxtrim = fractal_trimmer(vecN, cut)
  vecrtrim = vecr[idxtrim]
  vecNtrim = vecN[idxtrim]
  nnn      = length(vecrtrim)
  estdim   = sum(coef(lm(log(vecNtrim)~log(1/vecrtrim)))[2]) # lm fitting


  #   5-3. show and return
  output = list()
  output$estdim = estdim
  output$r  = vecr
  output$Nr = vecN
  return(output)
}

#' @keywords internal
#' @noRd
approx_boxcount <- function(r, Nr, nlevel){
  x = log(1/r)
  y = log(Nr)

  interp = stats::approx(x,y,n=nlevel)
  xnew   = interp$x
  ynew   = interp$y

  output = list()
  output$r = exp(-xnew)
  output$Nr = exp(ynew)
  return(output)
}


#  adjust log(counts) and returns index using percentage argument
#  returns the index
#' @keywords internal
#' @noRd
fractal_trimmer <- function(counts, cut){
  cut = sort(cut)
  counts = log(counts)
  fmin   = min(counts[is.finite(counts)])
  if (fmin <= 0){
    data = counts + abs(fmin)
  } else {
    data = counts
  }

  finite   = which(is.finite(data))
  maxvalue = max(data[is.finite(data)])
  thr1     = maxvalue*min(cut)
  thr2     = maxvalue*max(cut)
  idx = intersect(intersect(which(data>=thr1), which(data<=thr2)), finite)
  return(idx)
}
