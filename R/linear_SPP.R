#' Sparsity Preserving Projection
#'
#' Sparsity Preserving Projection (SPP) is an unsupervised linear dimension reduction technique.
#' It aims to preserve high-dimensional structure in a sparse manner to find projections
#' that keeps such sparsely-connected pattern in the low-dimensional space. Note that
#' we used \pkg{CVXR} for convenient computation, which may lead to slower execution
#' once used for large dataset.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param reltol tolerance level for stable computation of sparse reconstruction weights.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.integer(iris$Species)
#'
#' ## test different tolerance levels
#' out1 <- do.spp(X,ndim=2,reltol=0.001)
#' out2 <- do.spp(X,ndim=2,reltol=0.01)
#' out3 <- do.spp(X,ndim=2,reltol=0.1)
#'
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, main="SPP::reltol=.001")
#' plot(out2$Y, col=label, main="SPP::reltol=.01")
#' plot(out3$Y, col=label, main="SPP::reltol=.1")
#' par(opar)
#' }
#' @references
#' \insertRef{qiao_sparsity_2010}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_SPP
#' @concept linear_methods 
#' @export
do.spp <- function(X, ndim=2, preprocess = c("center","scale","cscale","decorrelate","whiten"), reltol=1e-4){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.spp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. reltol
  reltol = as.double(reltol)
  if (!check_NumMM(reltol,0,1,compact=FALSE)){stop("* do.spp : 'reltol' should be a small positive real number.")}

  #------------------------------------------------------------------------
  ## COMPUTATION PART 1 : PREPROCESSING
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION PART 2 : Main SPP : use CVXR
  #   1. compute S
  S = array(0,c(n,n))
  for (i in 1:n){
    #   1-1. separate vector and target matrix
    xi = matrix(pX[i,],nrow=p)
    Xi = t(pX[-i,])
    #   1-2. si index
    idxsi = rep(TRUE,n); idxsi[i] = FALSE
    #   1-3. compute si and return
    S[i,idxsi]=spp_compute_si(xi,Xi,reltol)
  }
  #   2. Sbeta
  Sbeta = S+t(S)-(t(S)%*%S)
  #   3. geigen with largest entries
  LHS = t(pX)%*%Sbeta%*%pX
  RHS = t(pX)%*%pX
  #   4. projection matrix : top vectors
  projection = aux.geigen(LHS,RHS,ndim,maximal=TRUE)

  #------------------------------------------------------------------------
  ## RETURN OUTPUT
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}

#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
spp_compute_si <- function(xi, Xi, reltol){
  n   = ncol(Xi)
  si  = CVXR::Variable(n)
  obj = CVXR::norm1(si)
  constr = list(CVXR::p_norm(xi-(Xi%*%si),p=2)<=reltol,(sum(si)==1))
  prob = CVXR::Problem(Minimize(obj), constr)
  result = solve(prob)
  return(as.vector(result$getValue(si)))
}
