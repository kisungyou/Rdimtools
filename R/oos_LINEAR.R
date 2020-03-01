#' Out-Of-Sample Prediction for Linear Methods
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
#' set.seed(46556)
#' idxtest  = sample(1:500,20)        # 20% of data for testing
#' idxtrain = setdiff(1:500,idxtest)  # 80% of data for training
#'
#' Xtrain = X[idxtrain,]
#' Xtest  = X[idxtest,]
#'
#' ## run PCA for train data
#' res_train = do.pca(Xtrain,ndim=2,preprocess="whiten")
#'
#' ## perform OOS.LINEAR on new dataset
#' ## note that inputs should be from a given model you trained
#' model.projection = res_train$projection
#' model.trfinfo    = res_train$trfinfo
#' res_test  = oos.linear(Xtest, model.projection, model.trfinfo)
#'
#' ## let's compare via visualization
#' xx = c(-2,2) # range of axis 1 for compact visualization
#' yy = c(-2,2) # range of axis 2 for compact visualization
#' mm = "black=train / red=test data" # figure title
#' YY = res_test$Ynew  # out-of-sample projection for test data
#'
#' plot(res_train$Y, type="p", xlim=xx, ylim=yy,
#'      main=mm, xlab="axis 1", ylab="axis 2")
#' points(YY[,1], YY[,2], lwd=3, col="red")
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
  if ((!('algtype'%in%names(trfinfo)))||(!('type'%in%names(trfinfo)))||(!('mean'%in%names(trfinfo)))||(!('multiplier'%in%names(trfinfo)))){
    stop("* oos.linear : 'trfinfo' is an invalid one.")
  }
  if (trfinfo$algtype!='linear'){
    stop("* oos.linear : this function is designed for 'linear' methods only.")
  }

  #------------------------------------------------------------------------
  # COMPUTE
  # 1. oos preprocessing
  Xtmp = aux.oospreprocess(Xnew, trfinfo)
  # 2. projection
  Y = Xtmp%*%projection

  #------------------------------------------------------------------------
  ## RETURN !
  result = list()
  result$Ynew = Y;
  return(result)
}



