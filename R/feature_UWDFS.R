#' Uncorrelated Worst-Case Discriminative Feature Selection
#'
#' Built upon \code{do.wdfs}, this method selects features step-by-step to opt out the redundant sets
#' by iteratively update feature scores via scaling by the correlation between target and previously chosen variables.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
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
#' set.seed(100)
#' subid    = sample(1:150,50)
#' iris.dat = as.matrix(iris[subid,1:4])
#' iris.lab = as.factor(iris[subid,5])
#'
#' ## compare with other algorithms
#' out1 = do.lda(iris.dat, iris.lab)
#' out2 = do.wdfs(iris.dat, iris.lab)
#' out3 = do.uwdfs(iris.dat, iris.lab)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=iris.lab, main="LDA")
#' plot(out2$Y, pch=19, col=iris.lab, main="WDFS")
#' plot(out3$Y, pch=19, col=iris.lab, main="UWDFS")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{liao_worstcase_2019}{Rdimtools}
#'
#'
#' @seealso \code{\link{do.wdfs}}
#' @rdname feature_UWDFS
#' @author Kisung You
#' @concept feature_methods
#' @export
do.uwdfs <- function(X, label, ndim=2, preprocess=c("null","center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label vector
  label  = check_label(label, n)
  ulabel = unique(label)
  C      = length(ulabel)
  if (C==1){
    stop("* do.uwdfs : 'label' should have at least 2 unique labelings.")
  }
  if (C==n){
    stop("* do.uwdfs : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.uwdfs : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PREPROCESSING OF THE DATA
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION : main part from WDFS
  #   1. class-mean and within-class scatter
  prepBC <- array(0,c(C,p))
  prepWC <- array(0,c(p,p,C))
  for (i in 1:C){
    idi = which(label==i)
    idn = length(idi)
    if (idn>1){
      prepBC[i,] = as.vector(base::colMeans(pX[idi,]))
    } else {
      prepBC[i,] = as.vector(pX[idi,])
    }
    prepWC[,,i] = stats::cov(pX[idi,])*(idn-1)/idn
  }
  #   2. compute scores
  wscore = rep(0,p)
  for (i in 1:p){
    wrvec = rep(0,p); wrvec[i] = 1
    wscore[i] = WDFS.score(wrvec, prepBC, prepWC)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : iterative update
  id.chosen = c()
  id.left   = (1:p)
  vec.score = wscore

  # main iteration
  for (i in 1:ndim){
    # select the one with the largest
    posmax = which.max(vec.score)
    idmax  = id.left[posmax]
    fmax   = as.vector(pX[,idmax])

    # update indices
    id.chosen = c(id.chosen, idmax)
    id.left   = id.left[-posmax]
    vec.score = vec.score[-posmax]

    # update scores : we only use Pearson
    for (j in 1:length(vec.score)){
      fj = pX[,id.left[j]]
      vec.score[j] = (1-base::abs(stats::cor(x=fmax, y=fj)))*vec.score[j]
    }
  }

  # convert to conventional format and compute projection
  idxvec = id.chosen
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$featidx = idxvec
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
