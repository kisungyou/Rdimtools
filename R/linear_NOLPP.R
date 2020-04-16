#' Nonnegative Orthogonal Locality Preserving Projection
#'
#' Nonnegative Orthogonal Locality Preserving Projection (NOLPP) is a variant of OLPP where
#' projection vectors - or, basis for learned subspace - contain no negative values.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param t kernel bandwidth in \eqn{(0,\infty)}.
#' @param maxiter number of maximum iteraions allowed.
#' @param reltol stopping criterion for incremental relative error.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \dontrun{
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])+50
#' label = as.integer(iris$Species)
#'
#' ## use different kernel bandwidths with 20% connectivity
#' out1 = do.nolpp(X, type=c("proportion",0.5), t=0.01)
#' out2 = do.nolpp(X, type=c("proportion",0.5), t=0.1)
#' out3 = do.nolpp(X, type=c("proportion",0.5), t=1)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, main="NOLPP::t=0.01")
#' plot(out2$Y, col=label, main="NOLPP::t=0.1")
#' plot(out3$Y, col=label, main="NOLPP::t=1")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zafeiriou_nonnegative_2010}{Rdimtools}
#'
#' @seealso \code{\link{do.olpp}}
#' @rdname linear_NOLPP
#' @author Kisung You
#' @concept linear_methods 
#' @export
do.nolpp <- function(X, ndim=2, type=c("proportion",0.1),
                     preprocess=c("null","center","scale","cscale","decorrelate","whiten"),
                     t=1.0, maxiter=1000, reltol=1e-5){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.nolpp : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   3. type
  nbdtype = type
  nbdsymmetric = "union"
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. t = kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t, .Machine$double.eps, Inf, compact=TRUE)){stop("* do.nolpp : 't' should be a bandwidth parameter in (0,Inf).")}
  #   * maxiter and reltol
  maxiter = as.integer(maxiter)
  if (!check_NumMM(maxiter, 5, 1e+6)){stop("* do.nolpp : 'maxiter' is a large positive integer for the number of iterations.")}
  reltol = as.double(reltol)
  if (!check_NumMM(reltol, .Machine$double.eps, 1)){stop("* do.nolpp : 'reltol' is a small positive real number for stopping criterion.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #   3. Dsqmat with kernelization
  Dsqmat = exp(-(as.matrix(dist(pX))^2)/t)

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR NONNEGATIVE OLPP
  #   1. compute auxiliary matrices
  A = Dsqmat*nbdmask
  L = diag(rowSums(A))-A

  #   2. compute cost function
  C = t(pX)%*%L%*%pX

  #   3. set initial matrix
  Uinit = matrix(runif(p*ndim),nrow=p)

  #   4. solve the minimization problem
  projection = method_nnprojmin(C, Uinit, reltol, maxiter)
  projection[(is.na(projection)||(is.infinite(projection)))] = 1
  for (i in 1:ndim){
    tgt = as.vector(projection[,i])
    projection[,i] = tgt/sqrt(sum(tgt^2))
  }
  projection = aux.adjprojection(projection)


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
