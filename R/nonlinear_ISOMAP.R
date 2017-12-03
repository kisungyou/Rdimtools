#' Isometric Feature Mapping
#'
#' \code{do.isomap} is an efficient implementation of a well-known \emph{Isomap} method
#' by Tenenbaum et al (2000). Its novelty comes from applying classical multidimensional
#' scaling on nonlinear manifold, which is approximated as a graph.
#'
#' @param X an \code{(n-by-p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.
#' @param weight \code{TRUE} to perform Isomap on weighted graph, or \code{FALSE} otherwise.
#' @param preprocess an additional option for preprocessing the data.
#' Default is ``center'' and other methods of ``decorrelate'', or ``whiten''
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \code{(n-by-ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#'@examples
#'## generate data
#'## in order to pass CRAN pretest, n is set to be small.
#'X <- aux.gensamples(n=28)
#'
#'## 1. connecting 10% of data for graph construction.
#'output1 <- do.isomap(X,ndim=2,type=c("proportion",0.10))
#'
#'## 2. constructing 25%-connected graph
#'output2 <- do.isomap(X,ndim=2,type=c("proportion",0.25))
#'
#'## 3. constructing 25%-connected with binarization
#'output3 <- do.isomap(X,ndim=2,type=c("proportion",0.25),weight=FALSE)
#'
#'## Visualize three different projections
#'par(mfrow=c(1,3))
#'plot(output1$Y[,1],output1$Y[,2],main="5%")
#'plot(output2$Y[,1],output2$Y[,2],main="10%")
#'plot(output3$Y[,1],output3$Y[,2],main="10%+Binary")
#'
#'
#'
#' @references
#' \insertRef{silva_global_2003}{Rdimtools}
#'
#' @rdname nonlinear_ISOMAP
#' @author Kisung You
#' @export
do.isomap <- function(X,ndim=2,type=c("proportion",0.1),symmetric="union",weight=TRUE,preprocess="center"){
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
  nbdsymmetric = symmetric
  if (!is.element(nbdsymmetric,c("union","intersect","asymmetric"))){
    stop("* do.isomap : 'symmetric' should have one of three values.")
  }
  algweight = weight
  if (!is.logical(algweight)){
    stop("* do.isomap : 'weight' should be a logical value.")
  }
  algpreprocess = preprocess
  if (!is.element(algpreprocess,c("center","whiten","decorrelate"))){
    stop("* do.isomap : 'preprocess' should be one of three values.")
  }

  # 3. process : data preprocessing
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  trfinfo$algtype = "nonlinear"
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
