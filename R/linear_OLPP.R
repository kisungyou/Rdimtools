#' Orthogonal Locality Preserving Projection
#'
#' Orthogonal Locality Preserving Projection (OLPP) is a variant of \code{do.lpp}, which
#' extracts orthogonal basis functions to reconstruct the data in a more intuitive fashion.
#' It adopts PCA as preprocessing step and uses only one eigenvector at each iteration in that
#' it might incur warning messages for solving near-singular system of linear equations.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}.
#' See also \code{\link{aux.graphnbd}} for more details.
#' @param weight \code{TRUE} to perform LPP on weighted graph, or \code{FALSE} otherwise.
#' @param preprocess an additional option for preprocessing the data. See \code{\link{aux.preprocess}} for details.
#' @param t bandwidth for heat kernel in \eqn{(0,\infty)}
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.integer(iris$Species)
#'
#' ##  connecting 10% and 25% of data for graph construction each.
#' output1 <- do.olpp(X,ndim=2,type=c("proportion",0.10))
#' output2 <- do.olpp(X,ndim=2,type=c("proportion",0.25))
#'
#' ## Visualize
#' #  In theory, it should show two separated groups of data
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(output1$Y, col=label, main="OLPP::10% connected")
#' plot(output2$Y, col=label, main="OLPP::25% connected")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{cai_orthogonal_2006}{Rdimtools}
#'
#' @seealso \code{\link{do.lpp}}
#' @author Kisung You
#' @rdname linear_OLPP
#' @export
do.olpp <- function(X,ndim=2,type=c("proportion",0.1),symmetric=c("union","intersect","asymmetric"),
                    weight=TRUE,preprocess=c("center","scale","cscale","decorrelate","whiten"),t=1.0){
  # Preprocessing : typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("*do.olpp : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # Preprocessing 2 : parameters
  # 2-1. aux.graphnbd
  #   type : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #   symmetric : 'intersect','union', or 'asymmetric'
  # 2-2. OLPP itself
  #   weight     : TRUE
  #   preprocess : 'null','center','decorrelate', or 'whiten'
  #   t          : heat kernel bandwidth

  nbdtype = type;
  nbdsymmetric = match.arg(symmetric)
  algweight = weight
  if (!is.logical(algweight)){
    stop("* do.olpp : 'weight' should be a logical variable.")
  }
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  if (!is.numeric(t)||(t<=0)||is.na(t)||is.infinite(t)){
    stop("* do.olpp : 't' is a bandwidth parameter in (0,infinity).")
  }

  # Preprocessing 3 : data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  ## MAIN COMPUTATION
  #   step 1. PCA preprocessing
  covX    = stats::cov(pX)
  eigtest = eigen(covX, only.values=TRUE)
  pcadim  = sum(eigtest$values > 0)

  if (pcadim <= ndim){
    warning("* do.olpp : target 'ndim' is larger than intrinsic data dimension achieved from PCA.")
    output = do.pca(X, ndim=pcadim,preprocess=algpreprocess)
    return(output)
  }

  pXpca = do.pca(pX, ndim=pcadim, preprocess="center")
  Xpca  = pXpca$Y
  Wpca  = aux.adjprojection(pXpca$projection)

  #   step 2. adjacency graph
  #   here, wD is now Distance Matrix, which is denoted as S in the note.
  nbdstruct = aux.graphnbd(Xpca, method="euclidean", type=nbdtype, symmetric=nbdsymmetric)
  D     = nbdstruct$dist
  Dmask = nbdstruct$mask
  nD    = ncol(D)
  # 5. process : nbd binarization
  wD = Dmask*D
  idnan = is.na(wD)
  wD[idnan] = 0
  if (!algweight){
    wD = wD*exp(-matrix(as.double(Dmask),nrow=nD)/t)
  }

  #   step 3. main projection matrix
  Wolpp = aux.adjprojection(method_olpp(t(Xpca), wD, ndim));

  #   step 4. computation !
  #   1. adjust projection matrix
  proj_all = (Wpca%*%Wolpp)

  #   2. return output
  result = list()
  result$Y = pX%*%proj_all
  result$projection = proj_all
  trfinfo$algtype = "linear"
  result$trfinfo = trfinfo
  return(result)
}
