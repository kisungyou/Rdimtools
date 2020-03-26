#' Conformal Isometric Feature Mapping
#'
#' Conformal Isomap(C-Isomap) is a variant of a celebrated method of Isomap. It aims at, rather than
#' preserving full isometry, maintaining infinitestimal angles - conformality - in that
#' it alters geodesic distance to reflect scale information.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.
#' @param weight \code{TRUE} to perform Isomap on weighted graph, or \code{FALSE} otherwise.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#'
#' @examples
#' \donttest{
#' ## generate data
#' X <- aux.gensamples(dname="cswiss",n=500)
#'
#' ## 1. original Isomap
#' output1 <- do.isomap(X,ndim=2)
#'
#' ## 2. C-Isomap
#' output2 <- do.cisomap(X,ndim=2)
#'
#' ## 3. C-Isomap on a binarized graph
#' output3 <- do.cisomap(X,ndim=2,weight=FALSE)
#'
#' ## Visualize three different projections
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(output1$Y, main="Isomap")
#' plot(output2$Y, main="C-Isomap")
#' plot(output3$Y, main="Binarized C-Isomap")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{silva_global_2003}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_CISOMAP
#' @export
do.cisomap <- function(X,ndim=2,type=c("proportion",0.1),symmetric=c("union","intersect","asymmetric"),weight=TRUE,
                       preprocess=c("center","scale","cscale","whiten","decorrelate")){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.cisomap : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # 2. ... parameters
  # 2-1. aux.graphnbd
  #   type : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #   symmetric : 'intersect','union', or 'asymmetric'
  # 2-2. isomap itself
  #   weight     : TRUE
  #   preprocess : 'center','decorrelate', or 'whiten'
  nbdtype = type
  if (missing(symmetric)){
    nbdsymmetric = "union"
  } else {
    nbdsymmetric = match.arg(symmetric)
  }


  algweight = weight
  if (!is.logical(algweight)){
    stop("* do.cisomap : 'weight' should be a logical variable.")
  }
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  # 3. process : data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  # 4. process : neighborhood selection
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  D     = nbdstruct$dist
  Dmask = nbdstruct$mask
  nD    = ncol(D)
  # 5-1. process : nbd binarization
  if (algweight){
    wD = Dmask*D
    idnan = is.na(wD)
    wD[idnan] = 0
  } else {
    wD = matrix(as.double(Dmask),nrow=nD)
  }

  # 5-2. conformality
  MDists = array(0,c(1,nD))
  for (i in 1:nD){
    MDists[i] = mean(D[i,Dmask[i,]])
  }
  scalermat   = diag(as.vector(1/sqrt(MDists)))
  wDconformal = (scalermat %*% wD %*% scalermat)

  # 6. process : shortest path
  sD = aux.shortestpath(wDconformal)

  # 7. main computation
  output = method_mdsD(sD);
  eigvals = rev(output$eigval)
  eigvecs = output$eigvec[,rev(seq_len(length(eigvals)))]

  # 8. output
  matS = diag(sqrt(eigvals[1:ndim]))
  matU = eigvecs[,1:ndim]
  result = list()
  result$Y = t(matS %*% t(matU));
  result$trfinfo = trfinfo
  return(result)
}
