#' Unsupervised Spectral Feature Selection
#'
#' SPEC algorithm selects features from the data via spectral graph approach.
#' Three types of ranking methods that appeared in the paper are available where
#' the graph laplacian is built via RBF kernel.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param sigma bandwidth parameter for RBK kernel of type \eqn{S_{i,j} = \exp(-\|x_i - x_j \|^2 / 2\sigma^2 )}.
#' @param ranking types of feature scoring method. See the paper in the reference for more details.
#' @param preprocess an additional option for preprocessing the data. Default is "null". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{sscore}{a length-\eqn{p} vector of spectral feature scores.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' ## it is known that feature 3 and 4 are more important.
#' data(iris)
#' iris.dat = as.matrix(iris[,1:4])
#' iris.lab = as.factor(iris[,5])
#'
#' ## try different ranking methods
#' mysig = 6
#' out1  = do.specu(iris.dat, sigma=mysig, ranking="method1")
#' out2  = do.specu(iris.dat, sigma=mysig, ranking="method2")
#' out3  = do.specu(iris.dat, sigma=mysig, ranking="method3")
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=iris.lab, main="SPECU::method1")
#' plot(out2$Y, col=iris.lab, main="SPECU::method2")
#' plot(out3$Y, col=iris.lab, main="SPECU::method3")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zhao_spectral_2007}{Rdimtools}
#'
#' @seealso \code{\link{do.specs}}
#' @rdname linear_SPECU
#' @author Kisung You
#' @concept linear_methods 
#' @export
do.specu <- function(X, ndim=2, sigma=1.0, ranking=c("method1","method2","method3"),
                     preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  # #   2. label vector
  # label  = check_label(label, n)
  #   3. ndim
  ndim = round(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.specu : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. other parameters
  myscore = match.arg(ranking)

  #------------------------------------------------------------------------
  ## COMPUTATION : Preliminary
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION : Graph Construction
  # index = list()
  # for (i in 1:length(unique(label))){
  #   index[[i]] = which(label==i)
  # }
  # indlength = rep(0,length(index))
  # for (i in 1:length(index)){
  #   indlength[i] = length(index[[i]])
  # }
  # affinity
  DM2 = (as.matrix(stats::dist(pX))^2)
  S   = exp(-DM2/(sigma^2))

  #------------------------------------------------------------------------
  ## COMPUTATION : normalized graph laplacian
  # normalized graph laplacian
  rowsumS  = base::rowSums(S)
  Dhalfvec = sqrt(rowsumS)
  Dhalfinv = diag(1/sqrt(rowsumS))
  L = diag(n) - Dhalfinv%*%S%*%Dhalfinv

  #------------------------------------------------------------------------
  ## COMPUTATION : case branching
  rankvec = switch(myscore,
                   "method1" = aux_spec_method1(pX, L, Dhalfvec),
                   "method2" = aux_spec_method2(pX, L, Dhalfvec),
                   "method3" = aux_spec_method3(pX, L, Dhalfvec))
  if (all(myscore=="method3")){
    idxvec = base::order(rankvec, decreasing=TRUE)[1:ndim]
  } else {
    idxvec = base::order(rankvec, decreasing=FALSE)[1:ndim]
  }
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$sscore  = rankvec
  result$featidx = idxvec
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
