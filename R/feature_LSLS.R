#' Locality Sensitive Laplacian Score
#'
#' Locality Sensitive Laplacian Score (LSLS) is a supervised linear feature extraction method that combines
#' a feature selection framework of laplacian score where the graph laplacian is adjusted as in the
#' scheme of LSDA. The adjustment is taken via decomposed affinity matrices which are separately constructed
#' using the provided class label information.
#'
#' @seealso \code{\link{do.lsda}}, \code{\link{do.lscore}}
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param alpha a weight factor; should be a real number in \eqn{[0,1]}.
#' @param k an integer; the size of a neighborhood.
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
#' ## compare different neighborhood sizes
#' out1 = do.lsls(iris.dat, iris.lab, k=3)
#' out2 = do.lsls(iris.dat, iris.lab, k=6)
#' out3 = do.lsls(iris.dat, iris.lab, k=9)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=iris.lab, pch=19, main="LSLS::k=3")
#' plot(out2$Y, col=iris.lab, pch=19, main="LSLS::k=6")
#' plot(out3$Y, col=iris.lab, pch=19, main="LSLS::k=9")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{liao_gene_2014a}{Rdimtools}
#'
#' @rdname feature_LSLS
#' @author Kisung You
#' @concept feature_methods
#' @export
do.lsls <- function(X, label, ndim=2, alpha=0.5, k=5, preprocess=c("null","center","scale","cscale","decorrelate","whiten")){
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
    stop("* do.lsls : 'label' should have at least 2 unique labelings.")
  }
  if (C==n){
    stop("* do.lsls : given 'label' has all unique elements.")
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.lsls : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. alpha (weight parameter in [0,1])
  myalpha = as.double(alpha)
  if ((length(alpha)>1)||(alpha<0)||(alpha>1)){
    stop("* do.lsls : weight factor 'alpha' should be a value in [0,1].")
  }
  #   6. k (neighborhood size)
  myk = round(k)

  #------------------------------------------------------------------------
  ## STEP 1. PREPROCESSING OF THE DATA
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  ## STEP 2. NEIGHBORHOOD INFORMATION
  nbdobj  = RANN::nn2(pX, k=(myk+1))
  nbdinfo = nbdobj$nn.idx[,2:(myk+1)]

  ## STEP 3. COMPUTE W, Wb, and Ww
  #  3-1. masking
  matW = method_lsls(pX, nbdinfo)
  vecD = base::rowSums(matW)
  #  3-2. separation
  matWb = array(0,c(n,n)) # between (different)
  matWw = array(0,c(n,n)) # within  (same)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (label[i]==label[j]){ # within (same)
        matWw[i,j] <- matWw[j,i] <- matW[i,j]
      } else {
        matWb[i,j] <- matWb[j,i] <- matW[i,j]
      }
    }
  }

  ## STEP 4. WEIGHTED LAPLACIAN (STILL MYSTERIOUS)
  matL = myalpha*(base::diag(base::rowSums(matWw))-matWw) - (1-myalpha)*(base::diag(base::rowSums(matWb))-matWb)

  ## STEP 5. COMPUTE LSLS SCORES
  lscore = rep(0,p)
  for (i in 1:p){
    vecr  = as.vector(pX[,i])                                      # choose the vector
    vecrc = vecr - (base::sum(vecr*vecD)/base::sum(vecD))*rep(1,n) # centering

    term1     = base::sum(vecrc*as.vector(matL%*%vecrc))
    term2     = base::sum(vecrc*(base::sum(vecrc*vecD)))
    lscore[i] = term1/term2
  }

  ## STEP 6. POST-PROCESSING
  #  6-1. choose the smallest ones
  idxvec = base::order(lscore, decreasing = FALSE)[1:ndim]
  #  6-2. projection matrix
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
