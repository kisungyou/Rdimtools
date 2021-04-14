#' Sammon Mapping
#'
#' \code{do.sammon} is an implementation for Sammon mapping, one of the earliest
#' dimension reduction techniques that aims to find low-dimensional embedding
#' that preserves pairwise distance structure in high-dimensional data space.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param initialize \code{"random"} or \code{"pca"}; the former performs
#' fast random projection (see also \code{\link{do.rndproj}}) and the latter
#' performs standard PCA (see also \code{\link{do.pca}}).
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @references Sammon, J.W. (1969) \emph{A Nonlinear Mapping for Data Structure Analysis}.
#' IEEE Transactions on Computers, C-18 5:401-409.
#'
#' @examples
#' \donttest{
#' ## load iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.factor(iris$Species)
#'
#' ## compare two initialization
#' out1 = do.sammon(X,ndim=2)                   # random projection
#' out2 = do.sammon(X,ndim=2,initialize="pca")  # pca as initialization
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(out1$Y, pch=19, col=label, main="out1:rndproj")
#' plot(out2$Y, pch=19, col=label, main="out2:pca")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{sammon_nonlinear_1969}{Rdimtools}
#'
#' @rdname nonlinear_SAMMON
#' @author Kisung You
#' @concept nonlinear_methods
#' @export
do.sammon <- function(X,ndim=2,preprocess=c("null","center","scale","cscale","decorrelate","whiten"),
                      initialize=c("pca","random")){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.sammon : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)

  # 2. ... parameters
  # 2-1. Sammon Mapping Itself
  #   preprocess : 'null', 'center','decorrelate', or 'whiten'
  #   initialize : 'random' (default as fast random projection) or 'pca'

  algpreprocess = match.arg(preprocess)
  initsammon = match.arg(initialize)

  # 3. preprocess
  #   3-1. centering or so
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   3-2. data reduction
  if (initsammon=="random"){
    tmpout = do.rndproj(pX,ndim=ndim,type='achlioptas')
    initY  = tmpout$Y
  } else {
    tmpout = do.pca(pX,ndim=ndim)
    initY = tmpout$Y
  }

  # 4. main computation
  dX = stats::dist(pX)
  if (any(as.vector(dX) <= 0)){
    result = list()
    result$Y = initY
    result$trfinfo  = trfinfo
  } else {
    result = list()
    result$Y = MASS::sammon(dX, y=initY, k=ndim, trace=FALSE)$points
    result$trfinfo  = trfinfo
  }
  return(result)
}
