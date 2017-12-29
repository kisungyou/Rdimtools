#' Locality Preserving Projection
#'
#' \code{do.lpp} is a linear approximation to Laplacian Eigenmaps. More precisely,
#' it aims at finding a linear approximation to the eigenfunctions of the Laplace-Beltrami
#' operator on the graph-approximated data manifold.
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
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center" and other options of "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param t bandwidth for heat kernel in \eqn{(0,\infty)}
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{eigval}{a vector of eigenvalues corresponding to basis expansion in an ascending order.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#'@examples
#'## generate data
#'## in order to pass CRAN pretest, n is set to be small.
#'X <- aux.gensamples(dname="twinpeaks",n=28)
#'
#'## 1. connecting 10% of data for graph construction.
#'output1 <- do.lpp(X,ndim=2,type=c("proportion",0.10))
#'
#'## 2. constructing 25%-connected graph
#'output2 <- do.lpp(X,ndim=2,type=c("proportion",0.25))
#'
#'## 3. constructing half-connected graph as binarized
#'output3 <- do.lpp(X,ndim=2,type=c("proportion",0.5),weight=FALSE)
#'
#'## Visualize three different projections
#'if ((!is.na(output1))&&(!is.na(output2))&&(!is.na(output3))){
#'par(mfrow=c(1,3))
#'plot(output1$Y[,1],output1$Y[,2],main="10%")
#'plot(output2$Y[,1],output2$Y[,2],main="25%")
#'plot(output3$Y[,1],output3$Y[,2],main="50%")
#'} else {
#'message("* do.lpp : example : one of three trials must have failed.")
#'}
#'
#'
#' @references
#' \insertRef{he_locality_2005}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_LPP
#' @export
do.lpp <- function(X,ndim=2,type=c("proportion",0.1),symmetric="union",weight=TRUE,preprocess="center",t=1.0){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("*do.lpp : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # 2. ... parameters
  # 2-1. aux.graphnbd
  #   type : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #   symmetric : 'intersect','union', or 'asymmetric'
  # 2-2. LPP itself
  #   weight     : TRUE
  #   preprocess : 'null','center','decorrelate', or 'whiten'
  #   t          : heat kernel bandwidth

  nbdtype = type;
  nbdsymmetric = symmetric
  if (!is.element(nbdsymmetric,c("union","intersect","asymmetric"))){
    stop("* do.lpp : flag 'symmetric' is invalid.")
  }
  algweight = weight
  if (!is.logical(algweight)){
    stop("* do.lpp : 'weight' should be a logical variable.")
  }
  algpreprocess = preprocess
  if (!is.element(algpreprocess,c("center","whiten","decorrelate"))){
    stop("* do.lpp : 'preprocess' should be one of 3 options.")
  }
  if (!is.numeric(t)||(t<=0)||is.na(t)||is.infinite(t)){
    stop("* do.lpp : 't' is a bandwidth parameter in (0,infinity).")
  }

  # 3. process : data preprocessing
  if (algpreprocess=="null"){
    trfinfo = list()
    trfinfo$type = "null"
    pX = as.matrix(X,nrow=nrow(X))
  } else {
    tmplist = aux.preprocess(X,type=algpreprocess)
    trfinfo = tmplist$info
    pX      = tmplist$pX
  }
  # 4. process : neighborhood selection
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
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

  # 6. main computation
  tpX = t(pX)
  output = tryCatch(
    {suppressWarnings(method_lpp(tpX,wD))},
    error=function(cond){
      message("* do.lpp : 'error' in the method detected.")
      return(NA)
    },
    warning=function(cond){
      message("* do.lpp : 'warning' in the method detected.")
      return(NA)
    }
  )

  eigvals = output$eigval
  eigvecs = output$eigvec
  if (dim(eigvecs)[1]==0){
    result = NA
    return(result)
  } else {
    # 7. return output
    #   1. adjust projection
    projection = aux.adjprojection(eigvecs[,1:ndim])
    #   2. return output
    result = list()
    result$Y = pX %*% projection
    result$eigval = eigvals[1:ndim]
    result$projection = projection
    trfinfo$algtype = "linear"
    result$trfinfo = trfinfo
    return(result)
  }
}
