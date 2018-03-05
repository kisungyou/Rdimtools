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
#' \dontrun{
#' ## generate default dataset
#' ## in order to pass CRAN pretest, n is set to be small.
#' X <- aux.gensamples(n=99)
#'
#' ## compare two initialization
#' out1 = do.sammon(X,ndim=2)                   # random projection
#' out2 = do.sammon(X,ndim=2,initialize="pca")  # pca as initialization
#'
#' par(mfrow=c(1,2))
#' plot(out1$Y[,1],out1$Y[,2],main="out1:rndproj")
#' plot(out2$Y[,1],out2$Y[,2],main="out2:pca")
#' }
#'
#' @references
#' \insertRef{sammon_nonlinear_1969}{Rdimtools}
#'
#'
#' @rdname nonlinear_SAMMON
#' @author Kisung You
#' @export
do.sammon <- function(X,ndim=2,preprocess=c("null","center","scale","cscale","decorrelate","whiten"),
                      initialize=c("random","pca")){
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
  tpX = t(pX)
  result = list()
  result$Y = method_sammon(tpX,initY)
  result$trfinfo  = trfinfo
  return(result)
}
