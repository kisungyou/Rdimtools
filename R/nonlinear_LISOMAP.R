#' Landmark Isometric Feature Mapping
#'
#' Landmark Isomap is a variant of Isomap in that
#' it first finds a low-dimensional embedding using a small portion of given dataset
#' and graft the others in a manner to preserve as much pairwise distance from
#' all the other data points to landmark points as possible.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param ltype on how to select landmark points, either \code{"random"} or \code{"MaxMin"}.
#' @param npoints the number of landmark points to be drawn.
#' @param preprocess an option for preprocessing the data. Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.
#' @param weight \code{TRUE} to perform Landmark Isomap on weighted graph, or \code{FALSE} otherwise.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \donttest{
#' ## generate data
#' ## in order to pass CRAN pretest, n is set to be small.
#' X <- aux.gensamples(n=123)
#'
#' ## 1. use 10/28 data points with quarter-connected graph
#' output1 <- do.lisomap(X,ndim=2,type=c("proportion",0.25),npoints=10)
#'
#' ## 2. use 10/28 data points with half-connected graph
#' output2 <- do.lisomap(X,ndim=2,type=c("proportion",0.50),npoints=10)
#'
#' ## 3. use 14/28 data points with half-connected graph
#' output3 <- do.lisomap(X,ndim=2,type=c("proportion",0.50),npoints=14)
#'
#' ## Visualize three different projections
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(output1$Y, main="5%")
#' plot(output2$Y, main="10%")
#' plot(output3$Y, main="50%")
#' par(opar)
#' }
#'
#' @seealso \code{\link{do.isomap}}
#' @references
#' \insertRef{silva_global_2003}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_LISOMAP
#' @concept nonlinear_methods 
#' @export
do.lisomap <- function(X,ndim=2,ltype=c("random","MaxMin"),npoints=max(nrow(X)/5,ndim+1),
                       preprocess=c("center","scale","cscale","decorrelate","whiten"),
                       type=c("proportion",0.1),symmetric=c("union","intersect","asymmetric"),weight=TRUE){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.lisomap : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # 2. ... parameters
  # 2-1. landmark selection
  #   ltype      : "random" (default) or "MaxMin"
  #   npoints    : (ndim+1 ~ nrow(X)/2)
  # 2-2. lmds itself
  #   preprocess : 'center','decorrelate', or 'whiten'
  # 2-3. aux.graphnbd
  #   type       : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #   symmetric  : 'intersect','union', or 'asymmetric'
  #   weight     : TRUE
  if (missing(ltype)){
    ltype = "random"
  } else {
    ltype = match.arg(ltype)
  }
  npoints = as.integer(round(npoints))
  if (!is.numeric(npoints)||(npoints<=ndim)||(npoints>(nrow(X)/2+1))||is.na(npoints)||is.infinite(npoints)){
    stop("* do.lisomap : the number of landmark points should be [ndim+1,#(total data points)/2].")
  }
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  nbdtype = type
  if (missing(symmetric)){
    nbdsymmetric = "union"
  } else {
    nbdsymmetric = match.arg(symmetric)
  }
  algweight = weight
  if (!is.logical(algweight)){
    stop("* do.lisomap : 'weight' param should be a logical value.")
  }

  # 3. Preprocess the data.
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

  # 7. select landmark points
  if (ltype=="random"){
    landmarkidx = sample(1:nrow(pX),npoints)
  } else if (ltype=="MaxMin"){
    landmarkidx = aux.MaxMinLandmark(sD,npoints,pdflag=TRUE)
  }
  if (length(landmarkidx)!=npoints){
    stop("* do.lisomap : landmark selection process is incomplete.")
  }

  # 5. MDS on landmark points
  sDlandmark = sD[landmarkidx,landmarkidx]
  output  = method_mdsD(sDlandmark);
  eigvals = rev(output$eigval)
  eigvecs = output$eigvec[,rev(seq_len(length(eigvals)))]

  idxnegs = which(eigvals<0)
  eigvals = eigvals[-idxnegs]
  eigvecs = eigvecs[,-idxnegs]

  tgtidx = max((as.integer(min(which((cumsum(eigvals)/sum(eigvals))>=0.8)))),2*ndim)
  matS   = diag(sqrt(eigvals[1:tgtidx]))
  matU   = eigvecs[,1:tgtidx]
  Lk     = (matS %*% t(matU));

  # 6. Distance-Based Triangulation
  #   6-1. pseudoinverse for mapping of (k-by-n) matrix Lk#
  Lksharp = array(0,c(nrow(Lk),ncol(Lk)))
  for (i in 1:nrow(Lk)){
    tgtvec = Lk[i,]
    lambda = sqrt(sum(tgtvec^2))
    Lksharp[i,] = tgtvec/(lambda^2)
  }
  #   6-2. pairwise distance matrix
  Deltan = (sD[landmarkidx,landmarkidx])^2
  deltamu = rowMeans(Deltan)
  #   6-3. Iterate over all data
  Ydbt = array(0,c(nrow(Lk),nrow(pX)))
  for (i in 1:nrow(pX)){
    deltax = (sD[i,landmarkidx])^2
    Ydbt[,i] = (Lksharp %*% (deltax-deltamu))/(-2)
  }

  # 7. PCA align
  tYdbt = t(Ydbt)
  pcaoutput = do.pca(tYdbt,ndim=ndim)

  # 8. return output
  result = list()
  result$Y = pcaoutput$Y
  result$trfinfo = trfinfo
  return(result)
}
