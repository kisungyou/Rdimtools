#' Constraint Score
#'
#' Constraint Score is a filter-type algorithm for feature selection using pairwise constraints.
#' It first marks all pairwise constraints as same- and different-cluster and
#' construct a feature score for both constraints. It takes ratio or difference of
#' feature score vectors and selects the indices with smallest values.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of class labels.
#' @param ndim an integer-valued target dimension.
#' @param score type of score measures from two score vectors of same- and different-class pairwise constraints; \code{"ratio"} and \code{"difference"} method. See the paper from the reference for more details.
#' @param lambda a penalty value for different-class pairwise constraints. Only valid for \code{"difference"} scoring method.
#' @param preprocess an additional option for preprocessing the data. Default is "null". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{cscore}{a length-\eqn{p} vector of constraint scores. Indices with smallest values are selected.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## use iris data
#' ## it is known that feature 3 and 4 are more important.
#' data(iris)
#' iris.dat = as.matrix(iris[,1:4])
#' iris.lab = as.factor(iris[,5])
#'
#' ## try different strategy
#' out1 = do.cscore(iris.dat, iris.lab, score="ratio")
#' out2 = do.cscore(iris.dat, iris.lab, score="difference", lambda=0)
#' out3 = do.cscore(iris.dat, iris.lab, score="difference", lambda=0.5)
#' out4 = do.cscore(iris.dat, iris.lab, score="difference", lambda=1)
#'
#' ## visualize
#' opar <- par(mfrow=c(2,2), no.readonly=TRUE)
#' plot(out1$Y, col=iris.lab, main="ratio")
#' plot(out2$Y, col=iris.lab, main="diff/lambda=0")
#' plot(out3$Y, col=iris.lab, main="diff/lambda=0.5")
#' plot(out4$Y, col=iris.lab, main="diff/lambda=1")
#' par(opar)
#'
#' @references
#' \insertRef{zhang_constraint_2008a}{Rdimtools}
#'
#' @seealso \code{\link{do.cscoreg}}
#' @rdname linear_CSCORE
#' @author Kisung You
#' @export
do.cscore <- function(X, label, ndim=2, score=c("ratio","difference"), lambda=0.5,
                      preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  m = nrow(X)
  n = ncol(X)
  #   2. label vector
  label  = check_label(label, m)
  #   3. ndim
  ndim = round(ndim)
  if (!check_ndim(ndim,n)){
    stop("* do.cscore : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. other parameters
  myscore = match.arg(score)
  mylbd   = as.double(lambda)

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN
  #   1. compute SM and SC matrix
  cmats = cscore_Smatrix(label)
  matSC = cmats$SC
  matSM = cmats$SM
  #   2. compute elementary vectors
  vecM = method_scoresum(pX, matSM)
  vecC = method_scoresum(pX, matSC)
  #   3. compute score according to the score type
  if (all(myscore=="ratio")){
    rankvec = vecM/vecC
  } else {
    rankvec = vecM - mylbd*vecC
  }
  #   4. select the smallest ones
  idxvec = base::order(rankvec, decreasing=FALSE)[1:ndim]
  #   5. find the projection matrix
  projection = aux.featureindicator(n,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$cscore  = rankvec
  result$featidx = idxvec
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}


# auxiliary function for CSCORE -------------------------------------------
#' @keywords internal
cscore_Smatrix <- function(label){
  n = length(label)
  matSM = array(0,c(n,n))
  matSC = array(0,c(n,n))
  for (i in 1:(n-1)){
    labi = label[i]
    for (j in (i+1):n){
      labj = label[j]
      if (labi==labj){
        matSM[i,j] <- matSM[j,i] <- 1
      } else {
        matSC[i,j] <- matSC[j,i] <- 1
      }
    }
  }
  mats = list()
  mats$SM = matSM
  mats$SC = matSC
  return(mats)
}
