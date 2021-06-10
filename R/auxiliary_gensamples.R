#' Generate model-based samples
#'
#' It generates samples from predefined shapes, set by \code{dname} parameter.
#' Also incorporated a functionality to add white noise with degree \code{noise}.
#'
#' @param n the number of points to be generated.
#' @param noise level of additive white noise.
#'
#' @param dname name of a predefined shape. Should be one of \describe{
#' \item{\code{"swiss"}}{swiss roll}
#' \item{\code{"crown"}}{crown}
#' \item{\code{"helix"}}{helix}
#' \item{\code{"saddle"}}{manifold near saddle point}
#' \item{\code{"ribbon"}}{ribbon}
#' \item{\code{"bswiss"}}{broken swiss}
#' \item{\code{"cswiss"}}{cut swiss}
#' \item{\code{"twinpeaks"}}{two peaks}
#' \item{\code{"sinusoid"}}{sinusoid on the circle}
#' \item{\code{"mobius"}}{mobius strip embedded in \eqn{\mathbf{R}^3}}
#' \item{\code{"R12in72"}}{12-dimensional manifold in \eqn{\mathbf{R}^{12}}}
#' }
#' @param ... extra parameters for the followings #' \tabular{lll}{
#' parameter \tab dname \tab description \cr
#' \code{ntwist} \tab \code{"mobius"} \tab number of twists
#' }
#'
#' @return an \eqn{(n\times p)} matrix of generated data by row. For all methods other than \code{"R12in72"}, it returns a matrix with \eqn{p=3}.
#'
#' @references
#' \insertRef{hein_intrinsic_2005}{Rdimtools}
#'
#' \insertRef{vandermaaten_learning_2009}{Rdimtools}
#'
#' @examples
#' \donttest{
#' ## generating toy example datasets
#' set.seed(100)
#' dat.swiss = aux.gensamples(50, dname="swiss")
#' dat.crown = aux.gensamples(50, dname="crown")
#' dat.helix = aux.gensamples(50, dname="helix")
#' }
#'
#' @author Kisung You
#' @rdname aux_gensamples
#' @export
aux.gensamples <- function(n=496,noise=0.01,
                           dname=c("swiss","crown","helix","saddle","ribbon","bswiss","cswiss","twinpeaks","sinusoid",
                                   "mobius","R12in72"), ...){
  #params = as.list(environment())
  n     = as.integer(n)
  noise = as.double(noise)
  dtype = (match.arg(dname))

  ## extra parameters
  extrapar = list(...)
  #   1. ntwist
  if ("ntwist" %in% names(extrapar)){
    ntwist = as.integer(extrapar$ntwist)
  } else {
    ntwist= as.integer(10)
  }
  if ((ntwist<1)||(is.na(ntwist))||(length(ntwist)>1)||(is.infinite(ntwist))){
    stop("* aux.gensamples : 'ntwist' should be a positive integer.")
  }

  # 2. type branching
  if (all(dtype=="swiss")){
    lowx = ((3*pi)/2)*(1+2*runif(n))

    highx = lowx*cos(lowx)
    highy = lowx*sin(lowx)
    highz = 30*runif(n)

    highD = cbind(highx,highy,highz) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (all(dtype=="crown")){
    lowD = (2*pi)*(1:n)/n

    highx = (2+sin(8*lowD))*sin(lowD)
    highy = (2+sin(8*lowD))*cos(lowD)
    highz = sin(8*lowD)

    highD = cbind(highx,highy,highz) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (all(dtype=="helix")){
    lowD = (2*pi)*(1:n)/n

    highx = (2+cos(8*lowD))*cos(lowD)
    highy = (2+cos(8*lowD))*sin(lowD)
    highz = sin(8*lowD)

    highD = cbind(highx,highy,highz) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (all(dtype=="saddle")){
    lowD = 1-2*matrix(rnorm(2*n),c(n,2))

    hx = lowD[,1]
    hy = lowD[,2]
    highz = sin(pi*hx)*tanh(3*hy)

    highD = cbind(lowD,highz) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (all(dtype=="ribbon")){
    lowD = (2*pi)*(1:n)/n

    hx = cos(lowD)
    hy = sin(lowD)
    h = 5*runif(n)

    highD = cbind(hx,hx*hy,h) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (all(dtype=="bswiss")){
    nceil = ceiling(n/2)
    lowD1 = (3*pi/2)*(1+0.8*runif(nceil))
    lowD2 = (2*pi/2)*(1+0.8*runif(n-nceil)+0.6)
    lowD  = c(lowD1,lowD2)

    height = 30*runif(n)

    hx = lowD*cos(lowD)
    hy = height
    hz = lowD*sin(lowD)

    highD = cbind(hx,hy,hz) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (all(dtype=="cswiss")){
    r = array(0,c(1,n))
    for (i in 1:n){
      pass = FALSE
      while (!pass){
        rr = runif(1)
        if (runif(1) > rr){
          r[i] = rr
          pass = TRUE
        }
      }
    }

    lowD = t((3*pi/2)*(1+(2*r)))
    height = 21*runif(n)

    hx = lowD*cos(lowD)
    hy = height
    hz = lowD*sin(lowD)

    highD = cbind(hx,hy,hz)+matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (all(dtype=="twinpeaks")){
    xi = runif(n)
    yi = runif(n)

    hx = 1-2*xi
    hy = sin(pi-2*pi*(xi))
    hz = tanh(3-6*yi)

    highD = cbind(hx,hy,hz)+matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (all(dtype=="sinusoid")){
    # hein audibert
    tt = runif(n, min=0, max=(2*3.141592))

    hx = sin(tt)
    hy = cos(tt)
    hz = sin(150*tt)/10

    highD = cbind(hx,hy,hz)+matrix(rnorm(n*3,sd=noise),nrow=n)
  } else if (all(dtype=="mobius")){
    u = runif(n,-1,1)
    v = runif(n,0,2*3.141592)

    hx = (1+(u/2)*cos(ntwist*v/2))*(cos(v))
    hy = (1+(u/2)*cos(ntwist*v/2))*(sin(v))
    hz = (u/2)*sin(ntwist*v/2)

    highD = cbind(hx,hy,hz)+matrix(rnorm(n*3,sd=noise),nrow=n)
  } else if (all(dtype=="R12in72")){
    xi = array(0,c(n,24)) # without noise
    for (it in 1:n){
      alpha = runif(12)
      for (i in 1:11){
        xi[it,(2*i-1)] = alpha[i+1]*cos(2*pi*alpha[i])
        xi[it,2*i]     = alpha[i+1]*sin(2*pi*alpha[i])
      }
      xi[it,23] = alpha[1]*cos(2*pi*alpha[12])
      xi[it,24] = alpha[1]*sin(2*pi*alpha[12])
    }

    highD = cbind(xi,xi,xi) + matrix(rnorm(n*72,sd=noise),nrow=n)
  }
  return(highD)
}
