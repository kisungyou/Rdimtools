#' Isometric Feature Mapping
#'
#' \code{do.isomap} is an efficient implementation of a well-known \emph{Isomap} method
#' by Tenenbaum et al (2000). Its novelty comes from applying classical multidimensional
#' scaling on nonlinear manifold, which is approximated as a graph.
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
#' @examples
#' ## generate data
#' set.seed(100)
#' X <- aux.gensamples(n=123)
#'
#' ## 1. connecting 10% of data for graph construction.
#' output1 <- do.isomap(X,ndim=2,type=c("proportion",0.10),weight=FALSE)
#'
#' ## 2. constructing 25%-connected graph
#' output2 <- do.isomap(X,ndim=2,type=c("proportion",0.25),weight=FALSE)
#'
#' ## 3. constructing 25%-connected with binarization
#' output3 <- do.isomap(X,ndim=2,type=c("proportion",0.50),weight=FALSE)
#'
#' ## Visualize three different projections
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(output1$Y, main="10%")
#' plot(output2$Y, main="25%")
#' plot(output3$Y, main="25%+Binary")
#' par(opar)
#'
#' @references
#' \insertRef{silva_global_2003}{Rdimtools}
#'
#' @rdname nonlinear_ISOMAP
#' @author Kisung You
#' @concept nonlinear_methods
#' @export
do.isomap <- function(X,ndim=2,type=c("proportion",0.1),symmetric=c("union","intersect","asymmetric"),
                      weight=FALSE,preprocess=c("center","scale","cscale","decorrelate","whiten")){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.isomap : 'ndim' is a positive integer in [1,#(covariates)].")
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
    stop("* do.isomap : 'weight' should be a logical value.")
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
  # 5. process : nbd binarization
  if (algweight){
    wD = Dmask*D
    idnan = is.na(wD)
    wD[idnan] = 0
  } else {
    wD = matrix(as.double(Dmask),nrow=nD)
  }
  # 6. process : shortest path
  sD = aux.shortestpath(wD)

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
