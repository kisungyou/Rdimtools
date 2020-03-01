#' Diffusion Maps
#'
#' \code{do.dm} discovers low-dimensional manifold structure embedded in high-dimensional
#' data space using Diffusion Maps (DM). It exploits diffusion process and distances in data space to find
#' equivalent representations in low-dimensional space.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param bandwidth a scaling parameter for diffusion kernel. Default is 1 and should be a nonnegative real number.
#' @param timescale a target scale whose value represents behavior of heat kernels at time \emph{t}. Default is 1 and should be a positive real number.
#' @param multiscale logical; \code{FALSE} is to use the fixed \code{timescale} value, \code{TRUE} to ignore the given value.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{eigvals}{a vector of eigenvalues for Markov transition matrix.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' ## generate swiss roll data
#' X <- aux.gensamples(n=200)
#'
#' ## compare different bandwidths
#' out1 <- do.dm(X,bandwidth=10)
#' out2 <- do.dm(X,bandwidth=100)
#' out3 <- do.dm(X,bandwidth=1000)
#'
#' ## visualize
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(out1$Y, main="DM::bandwidth=10")
#' plot(out2$Y, main="DM::bandwidth=100")
#' plot(out3$Y, main="DM::bandwidth=1000")
#' par(opar)
#' }
#'
#'@references
#'\insertRef{nadler_diffusion_2005}{Rdimtools}
#'
#'\insertRef{coifman_diffusion_2006}{Rdimtools}
#'
#' @rdname nonlinear_DM
#' @author Kisung You
#' @export
do.dm <- function(X,ndim=2,preprocess=c("null","center","scale","cscale","decorrelate","whiten"),
                  bandwidth=1.0,timescale=1.0,multiscale=FALSE){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  ndim = as.integer(ndim)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.dm : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  n = nrow(X)
  d = ncol(X)

  # 2. Parameters
  # 2-1. Common
  #   preprocess     : 'null'(default),'center','whiten','decorrelate'
  # 2-2. Diffusion Maps only
  #   bandwidth      : 1(default) or a real number >= 0
  #   timescale      : 1(default) or a real number >  0
  #   threshold      : 1e-7(default)

  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  if (!is.numeric(bandwidth)|(bandwidth<0)|is.infinite(bandwidth)){
    stop("* do.dm : 'bandwidth' should be a real number >= 0.")
  }
  if (!is.numeric(timescale)|(timescale<=0)|is.infinite(timescale)){
    stop("* do.dm : 'timescale' should be a positive real number > 0.")
  }
  if (!is.logical(multiscale)){
    stop("* do.dm : 'multiscale' should be a logical variable.")
  }
  if (multiscale==TRUE){
    message("* do.dm : when 'multiscale' is TRUE, the given timescale value is ignored.")
  }

  # 3. Preprocess
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  # 4. Main Computation : Scheme by Ann B. Lee (http://www.stat.cmu.edu/~annlee/software.htm)
  #   4-1. compute symmetric graph laplacian
  K = exp(-((as.matrix(dist(pX)))^2)/bandwidth)
  v = diag(1/sqrt(rowSums(K)))
  A = v%*%K%*%v
  #   4-2. SVD Computation
  svdA = base::svd(A)
  U    = svdA$u[,1:(ndim+1)]
  psi = U/matrix(rep(as.vector(U[,1]), (ndim+1)),ncol=(ndim+1),byrow=FALSE)
  phi = U*matrix(rep(as.vector(U[,1]), (ndim+1)),ncol=(ndim+1),byrow=FALSE)
  eigenvals = svdA$d

  #   4-3. compute embedding :: depending on timescale
  if (multiscale==FALSE){
    lambda_t = (svdA$d[2:(ndim+1)]^timescale)
    Y = psi[,2:(ndim+1)]*outer(rep(1,n),as.vector(lambda_t))
  } else {
    lambda_multi = eigenvals[2:(ndim+1)]/(1-eigenvals[2:(ndim+1)])
    Y = psi[,2:(ndim+1)]*outer(rep(1,n),as.vector(lambda_multi))
  }


  # 4. output
  result = list()
  result$Y = Y
  result$trfinfo = trfinfo
  result$eigvals = eigenvals[2:(ndim+1)]
  return(result)
}
