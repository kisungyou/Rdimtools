#' Out-of-sample prediction for linear methods
#'
#' Linear dimensionality reduction methods such as PCA, LPP, or ICA \emph{explicitly} returns a matrix
#' for mapping or projection. When we have new data, therefore, we can simply use the mapping provided.
#' Inputs \code{projection} and \code{trfinfo} should be brought from original model you trained.
#'
#' @param Xnew an \eqn{(m\times p)} matrix or data frame whose rows are observations. If a vector is given,
#' it will be considered as an \eqn{(1\times p)} matrix with single observation.
#' @param projection a \eqn{(p\times ndim)} projection matrix.
#' @param trfinfo a list containing transformation information generated from manifold learning algorithms.
#' See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Ynew}{an \eqn{(m\times ndim)} matrix whose rows are embedded observations.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate sample data and separate them
#' X = aux.gensamples(n=500)
#' idxtest  = sample(1:500,50)
#' idxtrain = setdiff(1:500,idxtest)
#'
#' Xtrain = X[idxtrain,]
#' Xtest  = X[idxtest,]
#'
#' ## run PCA for train data
#' res_train = do.pca(Xtrain,ndim=2,preprocess="whiten")
#'
#' ## perform OOS.LINEAR on new dataset
#' ## note that inputs should be from a given model you trained
#' res_test  = oos.linear(Xtest, res_train$projection, res_train$trfinfo)
#'
#' ## let's compare
#' par(mfrow=c(1,2))
#' plot(res_train$Y[,1], res_train$Y[,2], main="original PCA")
#' plot(res_test$Ynew[,1], res_test$Ynew[,2], main="OOS results")
#' }
#'
#' @author Kisung You
#' @rdname oos_LINEAR
#' @export
oos.linear <- function(Xnew, projection, trfinfo){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  # 0. if input Xnew is vector, make it a matrix.
  if (is.vector(Xnew)){
    Xnew = matrix(Xnew, nrow=1)
  }
  # 1. check Xnew and get size information
  aux.typecheck(Xnew)
  n = nrow(Xnew)
  p = ncol(Xnew)
  # 2. check projection
  aux.typecheck(projection)
  if (nrow(projection)!=p){
    stop("* oos.linear : 'projection' matrix must have ncol(Xnew) number of rows.")
  }
  ndim = ncol(projection)
  # 3. trfinfo
  if (!is.list(trfinfo)){
    stop("* oos.linear : 'trfinfo' should be provided a list.")
  }
  if ((!('algtype'%in%names(trfinfo)))||(!('type'%in%names(trfinfo)))){
    stop("* oos.linear : 'trfinfo' is an invalid one.")
  }
  if (trfinfo$algtype!='linear'){
    stop("* oos.linear : this function is designed for 'linear' methods only.")
  }

  #------------------------------------------------------------------------
  # COMPUTE
  Xtmp = array(0,c(n,p))
  type = trfinfo$type
  # 1. mean centering
  if (type!='null'){
    if (!('mean'%in%names(trfinfo))){
      stop("* oos.linear : 'trfinfo' should contain an entry called 'mean'.")
    }
    vecmean = as.vector(trfinfo$mean)
    if (length(vecmean)!=p){
      stop("* oos.linear : 'trfinfo$mean' has an invalid size.")
    }
    for (i in 1:n){
      Xtmp[i,]=Xnew[i,]-vecmean
    }
  }
  # 2. multiplication
  if ((type=="decorrelate")||(type=="whiten")){
    if (!('multiplier'%in%names(trfinfo))){
      stop("* oos.linear : 'trfinfo' should contain an entry called 'multiplier'.")
    }
    matmulti = trfinfo$multiplier
    aux.typecheck(matmulti)
    if ((nrow(matmulti)!=p)||(ncol(matmulti)!=p)){
      stop("* oos.linear : 'trfinfo$multiplier' has an invalid dimension.")
    }
    Xtmp = Xtmp%*%matmulti
  }
  # 3. projection
  Ynew = Xtmp%*%projection

  #------------------------------------------------------------------------
  ## RETURN !
  result = list()
  result$Ynew = Ynew;
  return(result)
}



