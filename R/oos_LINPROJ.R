#' OOS : Linear Projection
#'
#' The simplest way of out-of-sample extension might be linear regression even though the original embedding
#' is not the linear type by solving
#' \deqn{\textrm{min}_{\beta} \|Xold \beta - Yold\|_2^2} and use the estimate \eqn{\hat{beta}} to acquire
#' \deqn{Ynew = Xnew \hat{\beta}}. Due to the choice of original preprocessing, \code{trfinfo} must be brought
#' from the original model you trained.
#'
#' @param Xold an \eqn{(n\times p)} matrix of data in original high-dimensional space.
#' @param Yold an \eqn{(n\times ndim)} matrix of data in reduced-dimensional space.
#' @param trfinfo a list containing transformation information generated from manifold learning algorithms.
#' See also \code{\link{aux.preprocess}} for more details.
#' @param Xnew an \eqn{(m\times p)} matrix for out-of-sample extension.
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
#' Xold = X[idxtrain,]
#' Xnew = X[idxtest,]
#'
#' ## run PCA for train data
#' training = do.pca(Xold,ndim=2,preprocess="whiten")
#' Yold     = training$Y       # embedded data points
#' oldinfo  = training$trfinfo # preprocessing information
#'
#' ## perform out-of-sample extension
#' output  = oos.linproj(Xold, Yold, oldinfo, Xnew)
#' Ynew    = output$Ynew
#'
#' ## let's compare via visualization
#' xx = c(-2,2) # range of axis 1 for compact visualization
#' yy = c(-2,2) # range of axis 2 for compact visualization
#' mm = "black=train / red=test data" # figure title
#'
#' plot(Yold[,1], Yold[,2], type="p", xlim=xx, ylim=yy, main=mm, xlab="axis 1", ylab="axis 2")
#' par(new=TRUE)
#' plot(Ynew[,1], Ynew[,2], type="p", lwd=3, col="red", xlim=xx, ylim=yy, xlab="", ylab="")
#' }
#'
#' @author Kisung You
#' @rdname oos_LINPROJ
#' @export
oos.linproj <- function(Xold, Yold, trfinfo, Xnew){
  methodname = "linproj"
  #------------------------------------------------------------------------
  ## PREPROCESSING
  # 0. if inputs are vectors, make it a matrix.
  if (is.vector(Xnew)){    Xnew = matrix(Xnew, nrow=1)  }
  if (is.vector(Xold)){    Xold = matrix(Xold, nrow=1)  }
  if (is.vector(Yold)){    Yold = matrix(Yold, nrow=1)  }

  # 1. parameters
  n = nrow(Xold)
  p = ncol(Xold)
  ndim = ncol(Yold)
  m    = nrow(Xnew)

  if (nrow(Yold)!=n){
    stop(paste(" oos.",methodname," : 'Yold' should have same number of rows/observations as 'Xold'.",sep=""))
  }
  if (ndim >= p){
    stop(paste("* oos.",methodname," : non-sense alert : reduced dimensionality should've been smaller than original data dimension."))
  }
  if (ncol(Xnew)!=p){
    stop(paste("* oos.",methodname," : 'Xnew' should have same number of columns as 'Xold'.",sep=""))
  }

  # 2. trfinfo
  if (!is.list(trfinfo)){
    stop(paste("* oos.",methodname," : 'trfinfo' should be provided a list.",sep=""))
  }
  if ((!('algtype'%in%names(trfinfo)))||(!('type'%in%names(trfinfo)))||(!('mean'%in%names(trfinfo)))||(!('multiplier'%in%names(trfinfo)))){
    stop(paste("* oos.",methodname," : 'trfinfo' is an invalid one. Use 'info' output you acquired from one of functions in the package."))
  }


  #------------------------------------------------------------------------
  # COMPUTE : PREPROCESSING
  # old ones
  X = aux.oospreprocess(Xold, trfinfo)
  Y = Yold
  # new ones
  XX = aux.oospreprocess(Xnew, trfinfo)

  #------------------------------------------------------------------------
  # COMPUTE : MAIN PART FOR LINEAR PROJECTION
  # 1. X\beta = Y
  betahat = aux.pinv(X)%*%Y
  # 2. use the projection matrix for this
  YY      = XX%*%betahat


  #------------------------------------------------------------------------
  ## RETURN !
  result = list()
  result$Ynew = YY;
  return(result)
}


