#' Improved Local Tangent Space Alignment
#'
#' Conventional LTSA method relies on PCA for approximating local tangent spaces.
#' Improved LTSA (ILTSA) provides a remedy that can efficiently recover the geometric
#' structure of data manifolds even when data are sparse or non-uniformly distributed.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param t heat kernel bandwidth parameter in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \donttest{
#' ## generate ribbon-shaped data
#' X = aux.gensamples(dname="ribbon", n=200)
#'
#' ## try different bandwidth size
#' out1 <- do.iltsa(X, t=1)
#' out2 <- do.iltsa(X, t=100)
#' out3 <- do.iltsa(X, t=Inf)
#'
#' ## Visualize two comparisons
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="ILTSA::t=1")
#' plot(out2$Y, main="ILTSA::t=100")
#' plot(out3$Y, main="ILTSA::t=Inf")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zhang_improved_2011}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_ILTSA
#' @export
do.iltsa <- function(X, ndim=2, type=c("proportion",0.1),
                     symmetric=c("union","intersect","asymmetric"),
                     preprocess=c("center","scale","cscale","decorrelate","whiten"),
                     t=10.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.iltsa : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. nbd setup
  nbdtype = type
  if (missing(symmetric)){
    nbdsymmetric = "union"
  } else {
    nbdsymmetric = match.arg(symmetric)
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. t : kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t,.Machine$double.eps,Inf)){stop("* do.iltsa : 't' should be in (0,Inf).")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocess of data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #   2. build neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask


  #------------------------------------------------------------------------
  ## COMPUTATION : IMPROVED LTSA IN A DESCRIPTIVE ORDER
  #   1. main outer iterations
  C = array(0,c(n,n))
  for (i in 1:n){
    # 1-1. select the index
    Ii = which(nbdmask[i,])
    ki = length(Ii)
    if (ki <= 1){
      stop("* do.iltsa : select the larger neighborhood.")
    }
    # 1-2. compute Xi and wi
    tgtvec = as.vector(pX[i,])
    tgtmat = pX[Ii,]
    Xitilde = iltsa_Xtilde(tgtvec, tgtmat)
    vecWi   = iltsa_Wi(tgtvec, tgtmat, t)
    Wi      = diag(vecWi)
    # 1-3. compute orthonormal columns
    XWi    = (Xitilde%*%Wi)
    XWtXWi = (t(XWi)%*%XWi)
    eigXWi = RSpectra::eigs(XWtXWi, ndim)

    matEi = aux.adjprojection(eigXWi$vectors)
    vecdi = sqrt(eigXWi$values)
    # 1-4. Ti, Vi, and Ci
    Ti    = (diag(vecdi)%*%t(matEi)%*%diag(1/vecWi))
    invTi = aux.pinv(Ti)
    Vi    = (diag(ki)-(invTi%*%Ti))
    Ci    = (Vi%*%t(Vi))
    # 1-5. update scheme
    ek       = matrix(rep(1,ki),nrow=ki)
    C[Ii,Ii] = C[Ii,Ii]+Ci
    C[Ii,i]  = C[Ii,i]-as.vector(Ci%*%ek)
    C[i,Ii]  = C[i,Ii]-as.vector(t(ek)%*%Ci)
    C[i,i]   = C[i,i]+ as.double(t(ek)%*%Ci%*%ek)
  }

  #   2. compute (ndim+1) lowest eigenvectors
  if (ndim==(p-1)){
    Youtput = aux.adjprojection(base::eigen(C)$vectors[,p:2])
  } else {
    Youtput = aux.adjprojection(RSpectra::eigs(C, (ndim+1), which="SR")$vectors[,(ndim+1):2])
  }

  #------------------------------------------------------------------------
  ## RETURN OUTPUT
  result = list()
  result$Y = Youtput
  result$trfinfo = trfinfo
  return(result)
}





# . -----------------------------------------------------------------------
#' @keywords internal
#' @noRd
iltsa_Xtilde <- function(vec, mat){
  k = nrow(mat)
  n = ncol(mat)
  X = array(0,c(n,k))
  vec1 = as.vector(vec)
  for (i in 1:k){
    vec2  = as.vector(mat[i,])
    X[,i] = vec2-vec1
  }
  return(X)
}
#' @keywords internal
#' @noRd
iltsa_Wi <- function(vec, mat, t){
  k = nrow(mat)
  n = ncol(mat)
  output = rep(0,k)
  for (i in 1:k){
    vecdiff   = as.vector(vec)-as.vector(mat[i,])
    output[i] = exp(-sum(vecdiff*vecdiff)/t)
  }
  return(output)
}
