#' Landmark Multidimensional Scaling
#'
#' Landmark MDS is a variant of Classical Multidimensional Scaling in that
#' it first finds a low-dimensional embedding using a small portion of given dataset
#' and graft the others in a manner to preserve as much pairwise distance from
#' all the other data points to landmark points as possible.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param ltype on how to select landmark points, either \code{"random"} or \code{"MaxMin"}.
#' @param npoints the number of landmark points to be drawn.
#' @param preprocess an option for preprocessing the data. Default is "center".
#' See also \code{\link{aux.preprocess}} for more details.
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
#' # generate data
#' X <- aux.gensamples(dname="crown")
#'
#' ## 1. use 10% of random points
#' output1 <- do.lmds(X,ndim=2,npoints=round(nrow(X)/10))
#'
#' ## 2. using MaxMin scheme
#' output2 <- do.lmds(X,ndim=2,npoints=round(nrow(X)/10),ltype="MaxMin")
#'
#' ## 3. original mds case
#' output3 <- do.mds(X,ndim=2)
#'
#' ## Visualization
#' par(mfrow=c(1,3))
#' plot(output1$Y[,1],output2$Y[,2],main="10% random points")
#' plot(output2$Y[,1],output2$Y[,2],main="10% MaxMin points")
#' plot(output3$Y[,1],output3$Y[,2],main="original MDS")
#' }
#'
#' @seealso \code{\link{do.mds}}
#' @references
#' \insertRef{silva_global_2002}{Rdimtools}
#'
#' \insertRef{lee_landmark_2009}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LMDS
#' @export
do.lmds <- function(X,ndim=2,ltype="random",npoints=max(nrow(X)/5,ndim+1),
                    preprocess=c("center","cscale","decorrelate","whiten")){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.lmds : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # 2. ... parameters
  # 2-1. landmark selection
  #   ltype      : "random" (default) or "MaxMin"
  #   npoints   : (ndim+1 ~ nrow(X)/2)
  # 2-2. lmds itself
  #   preprocess : 'center','decorrelate', or 'whiten'
  if (!is.element(ltype,c("random","MaxMin"))){
    stop("* do.lmds : 'ltype' is either 'random' or 'MaxMin'.")
  }
  npoints = as.integer(round(npoints))
  if (!is.numeric(npoints)||(npoints<=ndim)||(npoints>nrow(X)/2)||is.na(npoints)||is.infinite(npoints)){
    stop("* do.lmds : the number of landmark points should be [ndim+1,#(total data points)/2].")
  }
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  # 3. Preprocess the data.
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX


  # 4. select landmark points
  if (ltype=="random"){
    landmarkidx = sample(1:nrow(pX),npoints)
  } else if (ltype=="MaxMin"){
    landmarkidx = aux.MaxMinLandmark(pX,npoints)
  }
  if (length(landmarkidx)!=npoints){
    stop("* do.lmds : landmark selection process is incomplete.")
  }

  # 5. MDS on landmark points
  pcarun = do.pca(pX[landmarkidx,])
  Lk = t(pcarun$Y)
  if (nrow(Lk)<=ndim){
    pcarun = do.pca(pX[landmarkidx,],ndim=ndim)
    Lk = t(pcarun$Y)
  }

  # 6. Distance-Based Triangulation
  pD = as.matrix(dist(pX))
  #   6-1. pseudoinverse for mapping of (k-by-n) matrix Lk#
  Lksharp = array(0,c(nrow(Lk),ncol(Lk)))
  for (i in 1:nrow(Lk)){
    tgtvec = Lk[i,]
    lambda = sqrt(sum(tgtvec^2))
    Lksharp[i,] = tgtvec/(lambda^2)
  }
  #   6-2. pairwise distance matrix
  Deltan = (pD[landmarkidx,landmarkidx])^2
  deltamu = rowMeans(Deltan)
  #   6-3. Iterate over all data
  Ydbt = array(0,c(nrow(Lk),nrow(pX)))
  for (i in 1:nrow(pX)){
    deltax = (pD[i,landmarkidx])^2
    Ydbt[,i] = (Lksharp %*% (deltax-deltamu))/(-2)
  }

  # 7. PCA align
  tYdbt = t(Ydbt)
  pcaoutput = do.pca(tYdbt,ndim=ndim,preprocess = "center")

  # 8. return output
  result = list()
  result$Y = pcaoutput$Y   # Y
  result$trfinfo = trfinfo # trfinfo

  LHS = t(pX)%*%pX
  RHS = t(pX)%*%(result$Y)
  result$projection = aux.adjprojection(solve(LHS,RHS)) # projection
  return(result)
}
