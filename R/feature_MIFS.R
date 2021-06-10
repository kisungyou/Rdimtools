#' Mutual Information for Selecting Features
#'
#' MIFS is a supervised feature selection that iteratively increases the subset of variables by choosing maximally informative feature based on the mutual information.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of class labels.
#' @param ndim an integer-valued target dimension.
#' @param beta penalty for relative importance of mutual information between the candidate and already-chosen features in iterations. Author proposes to use a value in \eqn{(0.5,1)}.
#' @param discretize the method for each variable to be discretized. The paper proposes \code{"default"} method to use 10 bins while \code{"histogram"} uses automatic discretization via Sturges' method.
#' @param preprocess an additional option for preprocessing the data. Default is "null". See also \code{\link{aux.preprocess}} for more details.
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
#' iris.dat = as.matrix(iris[,1:4])
#' iris.lab = as.factor(iris[,5])
#'
#' ## try different beta values
#' out1 = do.mifs(iris.dat, iris.lab, beta=0)
#' out2 = do.mifs(iris.dat, iris.lab, beta=0.5)
#' out3 = do.mifs(iris.dat, iris.lab, beta=1)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=iris.lab, main="beta=0")
#' plot(out2$Y, pch=19, col=iris.lab, main="beta=0.5")
#' plot(out3$Y, pch=19, col=iris.lab, main="beta=1")
#' par(opar)
#' }
#'
#'
#' @references
#' \insertRef{battiti_using_1994}{Rdimtools}
#'
#' @rdname feature_MIFS
#' @author Kisung You
#' @concept feature_methods
#' @export
do.mifs <- function(X, label, ndim=2, beta=0.75, discretize=c("default","histogram"),
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
    stop("* do.mifs : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. other parameters
  if (missing(discretize)){
    mydisc = "default"
  } else {
    mydisc = match.arg(discretize)
  }
  mybeta = as.double(beta)
  # if ((mybeta < 0.5)||(mybeta > 1)){
  #   message("* do.mifs : 'beta' is suggested to be in (0.5,1) from the author.")
  # }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #  Preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #  Discretization into labels
  if (all(mydisc=="default")){
    labX = base::apply(pX, 2, mifs_cut)
  } else {
    labX = base::apply(pX, 2, mifs_hist)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION
  #  Initial
  vecF = 1:n
  tmpF = rep(0,n)
  for (i in 1:n){
    tmpF[i] = as.double(mclustcomp::mclustcomp(label, labX[,i], types="mi")[2])
  }
  fmax = which.max(tmpF)
  vecF = base::setdiff(vecF, fmax)
  vecS = fmax

  #  Iterate
  MIcf = rep(0,n)
  MIfs = array(NA, c(n,n))
  for (it in 1:(ndim-1)){
    # how many candidates do we have ?
    nf = length(vecF)
    ns = (n-nf)
    # 1. compute I(c;f) if not available
    for (i in 1:nf){
      idf = vecF[i]
      if (is.na(MIcf[idf])){
        MIcf[idf] = as.double(mclustcomp::mclustcomp(label, labX[,idf],types="mi")[2])
      }
    }
    # 2. compute I(f;s) if not available
    for (i in 1:nf){
      idf = vecF[i]
      for (j in 1:ns){
        ids = vecS[j]
        if (is.na(MIfs[idf,ids])){
          MIfs[idf,ids] = as.double(mclustcomp::mclustcomp(labX[,idf],labX[,ids],types="mi")[2])
          MIfs[ids,idf] = MIfs[idf,ids]
        }
      }
    }
    # 3. now choose f that maximizes the cost function
    tmpcost = rep(0,nf)
    for (i in 1:nf){
      idf = vecF[i]
      term1 = MIcf[idf]
      term2 = (mybeta)*base::sum(MIfs[idf,vecS])
      tmpcost[i] = term1-term2
    }
    fmax = vecF[which.max(tmpcost)]
    # 4. update
    vecF = base::setdiff(vecF, fmax)
    vecS = c(vecS, fmax)
  }


  #  select the index and find the projection
  idxvec = vecS
  projection = aux.featureindicator(n,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$featidx = idxvec
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}


# auxiliary functions -----------------------------------------------------
#' @keywords internal
mifs_cut <- function(myvec){
  group = as.integer(base::cut(myvec, 10)) # Kf = 10 by default
  return(group)
}
#' @keywords internal
mifs_hist <- function(myvec){
  hobj = graphics::hist(myvec, plot=FALSE)
  return(base::findInterval(myvec, hobj$breaks))
}
