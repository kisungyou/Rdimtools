#' Collaborative Representation-based Projection
#'
#' Collaborative Representation-based Projection (CRP) is an unsupervised linear
#' dimension reduction method. Its embedding is based on \eqn{\ell}_2 graph construction,
#' similar to that of SPP where sparsity constraint is imposed via \eqn{\ell_1} optimization problem.
#' Note that though it may be way faster, rank deficiency can pose a great deal of problems,
#' especially when the dataset is large.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param lambda regularization parameter for constructing \eqn{\ell_2} graph.
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
#' ## generate samples
#' X <- aux.gensamples(n=200)
#'
#' ## test different regularization parameters
#' out1 <- do.crp(X,ndim=2,lambda=0.1)
#' out2 <- do.crp(X,ndim=2,lambda=1)
#' out3 <- do.crp(X,ndim=2,lambda=10)
#'
#' # visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="lambda=0.1")
#' plot(out2$Y[,1], out2$Y[,2], main="lambda=1")
#' plot(out3$Y[,1], out3$Y[,2], main="lambda=10")
#' }
#'
#' @seealso \code{\link{do.spp}}
#' @author Kisung You
#' @rdname linear_CRP
#' @export
do.crp <- function(X, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"), lambda=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.crp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. lambda
  lambda = as.double(lambda)
  if (!check_NumMM(lambda,0,1e+10,compact=TRUE)){stop("* do.crp : 'lambda' should be a nonnegative regularization parameter.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing the data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. PCA preprocessing
  eigtest = eigen(cov(pX), only.values=TRUE)
  pcadim  = sum(eigtest$values > 0)
  if (pcadim <= ndim){
    warning("* do.crp : target 'ndim' is larger than intrinsic data dimension achieved from PCA.")
    projection_first = diag(p)
    pcapX = pX
  } else{
    projection_first = aux.adjprojection(eigen(cov(pX))$vectors[,1:pcadim])
    pcapX = pX%*%projection_first
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR CRP
  #   1. compute Weight matrix
  W = array(0,c(n,n))
  seq1n = (1:n)
  for (i in 1:n){
    #   1-1. index
    currentidx = setdiff(seq1n,i)
    #   1-2. separate
    matX = t(pcapX[currentidx,])
    vecxi= as.vector(pcapX[i,])
    #   1-3. solve !
    solwi = crp_l2weight(matX, vecxi, lambda)
    #   1-4. assign
    W[i,currentidx] = solwi
  }
  #   2. build total scatter and local scatter
  allmean = colMeans(pcapX)
  St = aux_scatter(pcapX, allmean)
  Sl = t(pcapX)%*%(diag(n)-W-t(W)+(W%*%t(W)))%*%pcapX
  #   3. solve for matrix inversion
  solinv = aux.bicgstab(Sl, St, verbose=FALSE)$x
  #   4. extract projection eigenvectors
  projection_second = aux.adjprojection(RSpectra::eigs(solinv, ndim)$vectors)

  #------------------------------------------------------------------------
  ## RETURN RESULT
  #   1. combine two projection
  projection = projection_first%*%projection_second
  #   2. return
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}



#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
crp_l2weight <- function(A,y,lambda){
  # overdetermined and underdetermined case..
  #   case 1. underdetermined
  m = nrow(A)
  n = ncol(A)
  if (nrow(A)<ncol(A)){
    invinside = A%*%t(A)+lambda*diag(m)
    solreg    = as.vector(aux.bicgstab(invinside, y, verbose=FALSE)$x)
    sol       = as.vector(t(A)%*%matrix(solreg))
  } else {
    LHS = (t(A)%*%A + lambda*diag(n))
    RHS = as.vector(t(A)%*%y)
    sol = as.vector(aux.bicgstab(LHS, RHS, verbose=FALSE)$x)
  }
  return(sol)
}
