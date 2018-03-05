#' Locally Linear Embedding
#'
#' Locally-Linear Embedding (LLE) was introduced approximately at the same time as Isomap.
#' Its idea was motivated to describe entire data manifold by making a chain of local patches
#' in that low-dimensional embedding should resemble the connectivity pattern of patches.
#' \code{do.lle} also provides an automatic choice of regularization parameter based on an
#' optimality criterion suggested by authors.
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}.
#' See also \code{\link{aux.graphnbd}} for more details.
#' @param weight \code{TRUE} to perform LLE on weighted graph, or \code{FALSE} otherwise.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param regtype \code{TRUE} for automatic regularization parameter selection, \code{FALSE} otherwise as default.
#' @param regparam regularization parameter.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{eigvals}{a vector of eigenvalues from computation of embedding matrix.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Generate Data
#' X = aux.gensamples(n=100)
#'
#' ## 1. connecting 10% of data for graph construction.
#' output1 <- do.lle(X,ndim=2,type=c("proportion",0.10))
#'
#' ## 2. constructing 20%-connected graph
#' output2 <- do.lle(X,ndim=2,type=c("proportion",0.20))
#'
#' ## 3. constructing 50%-connected with bigger regularization parameter
#' output3 <- do.lle(X,ndim=2,type=c("proportion",0.5),regparam=10)
#'
#' ## Visualize three different projections
#' par(mfrow=c(1,3))
#' plot(output1$Y[,1],output1$Y[,2],main="5%")
#' plot(output2$Y[,1],output2$Y[,2],main="10%")
#' plot(output3$Y[,1],output3$Y[,2],main="50%+Binary")
#' }
#'
#' @seealso \href{https://www.cs.nyu.edu/~roweis/lle/}{Prof.Roweis' website}
#' @references
#' \insertRef{roweis_nonlinear_2000}{Rdimtools}
#'
#'
#' @author Kisung You
#' @rdname nonlinear_LLE
#' @export
do.lle <- function(X,ndim=2,type=c("proportion",0.1),symmetric="union",weight=TRUE,
                   preprocess=c("null","center","scale","cscale","decorrelate","whiten"),
                   regtype=FALSE, regparam=1.0){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("ERROR : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # 2. ... parameters
  # 2-1. aux.graphnbd
  #   type : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #   symmetric : 'intersect','union', or 'asymmetric'
  # 2-2. LLE itself
  #   weight     : TRUE
  #   preprocess : 'null', 'center','decorrelate', or 'whiten'
  #   regtype    : FALSE (default; need a value) / TRUE
  #   regparam   : 1 (default)

  nbdtype = type
  nbdsymmetric = symmetric
  if (!is.element(nbdsymmetric,c("intersect","union","asymmetric"))){
    stop("* do.lle : 'symmetric' option should have one of three types.")
  }

  algweight = TRUE
  if (!is.logical(algweight)){
    stop("* do.lle : 'weight' is a logical variable.")
  }
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  if (!is.logical(regtype)){
    stop("* do.lle : 'regtype' should be a logical variable.")
  }
  if (!is.numeric(regparam)||is.na(regparam)||is.infinite(regparam)||(regparam<=0)){
    stop("* do.lle : 'regparam' should be a positive real-valued number; it is a Tikhonov Factor.")
  }

  #   regtype    : FALSE (default; need a value) / TRUE
  #   regparam   : 1 (default)

  # 3. process : data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  n = nrow(pX)
  p = ncol(pX)

  # 4. process : neighborhood selection
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)


  # 5. main 1 : compute Weights
  #   k = max(apply(nbdstruct$mask,2,function(x) sum(as.double((x==TRUE)))))
  W = array(0,c(n,n))
  if (regtype==TRUE){
    regvals = array(0,c(1,n))
  }
  for (i in 1:n){
    #   5-1. separate target mask vector
    tgtmask = nbdstruct$mask[i,]
    tgtidx  = which(tgtmask==TRUE)

    #   5-2. select data
    #        For convenience, target matrix is transposed for Armadillo
    vec_tgt = pX[i,]
    mat_tgt = t(pX[tgtidx,])
    k = ncol(mat_tgt)

    #   5-3. compute with regularization
    #   5-3-1. No Automatic Regularization
    if (regtype==FALSE){
      w = method_lleW(mat_tgt,vec_tgt,regparam);
    } else {
      #   5-3-2. Automatic Regularization
      outW = method_lleWauto(mat_tgt,vec_tgt);
      w          = outW$w
      regvals[i] = outW$regparam
    }
    W[i,tgtidx] = w;
  }

  # 6. Main 2 : Compute Low-Dimensional Embedding
  embedding = method_lleM(W);

  # 7. Output
  #   this uses lowest (ndim+1) eigenpairs
  eigvals = embedding$eigval
  eigvecs = embedding$eigvec
  idxstart = max(min(which(eigvals>0)),2)

  result = list()
  result$Y = eigvecs[,idxstart:(idxstart+ndim-1)]
  result$trfinfo = trfinfo
  result$eigvals = eigvals[idxstart:(idxstart+ndim-1)]
  return(result)
}
