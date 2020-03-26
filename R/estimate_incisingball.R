#' Intrinsic Dimension Estimation with Incising Ball
#'
#' Incising ball methods exploits the exponential relationship of the number of samples
#' contained in a ball and the radius of the incising ball.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{estimated intrinsic dimension.}
#' }
#'
#' @examples
#' \donttest{
#' ## create an example data with intrinsic dimension 2
#' X = cbind(aux.gensamples(dname="swiss"),aux.gensamples(dname="swiss"))
#'
#' ## acquire an estimate for intrinsic dimension
#' output = est.incisingball(X)
#' sprintf("* est.incisingball : estimated dimension is %d.",output$estdim)
#' }
#'
#' @references
#' \insertRef{fan_intrinsic_2009}{Rdimtools}
#'
#' @rdname estimate_incisingball
#' @author Kisung You
#' @export
est.incisingball <- function(X){
  #####################################
  # Preprocessing
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  d = stats::dist(X)
  D = as.matrix(d)

  #####################################
  # STEP 1 : Find R and r
  dhist = hist(as.vector(d))
  ddensity = dhist$density
  dbreaks  = dhist$breaks

  didx = which(ddensity==max(ddensity))
  if (length(didx)>1){
    didx = didx[length(didx)] # find the largest one
  }
  R = (dbreaks[didx]+dbreaks[didx+1])/2
  r = min(as.vector(d))

  #####################################
  # STEP 2 : compute x and y
  x = rep(0,n)
  y = rep(0,n)
  for (i in 1:n){
    xi   = r + (i*(R-r)/n)
    x[i] = xi
    y[i] = (sum(D<xi))/n
  }

  #####################################
  # STEP 3 : least square fitting to find polynomials
  # I'll use C = 10^3 = 1000, alpha = 0.01 for testing beta_m
  order = min(n,p)
  oflag = FALSE
  while (oflag==FALSE){
    oflag = aux_inciseorder(y,x,order) # run it
    order = order - 1
    if (order==1){
      oflag=TRUE
    }
  }
  if (order==0){
    order = 1
  }

  #####################################
  # STEP 4 : return the result
  output = list()
  output$estdim = order
  return(output)
}


#' @keywords internal
#' @noRd
aux_inciseorder <- function(y,x,order){
  sumlm = summary(lm(y~poly(x,degree=order)))
  sumcoeff = matrix(sumlm$coefficients, ncol=4)

  bm = sumcoeff[order+1,1] # current
  bothers = sumcoeff[1:order,1]
  bmpval = sumcoeff[order+1,4]

  cond1 = (bm > 0)
  cond2 = (abs(max(bothers)/bm) < 1000)
  cond3 = (bmpval < 0.01)
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
