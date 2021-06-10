#' OOS : Linear Projection
#'
#' The simplest way of out-of-sample extension might be linear regression even though the original embedding
#' is not the linear type by solving
#' \deqn{\textrm{min}_{\beta} \|X_{old} \beta - Y_{old}\|_2^2} and use the estimate \eqn{\hat{beta}} to acquire
#' \deqn{Y_{new} = X_{new} \hat{\beta}}.
#'
#' @param Xold an \eqn{(n\times p)} matrix of data in original high-dimensional space.
#' @param Yold an \eqn{(n\times ndim)} matrix of data in reduced-dimensional space.
#' @param Xnew an \eqn{(m\times p)} matrix for out-of-sample extension.
#'
#' @return an \eqn{(m\times ndim)} matrix whose rows are embedded observations.
#'
#' @examples
#' \donttest{
#' ## generate sample data and separate them
#' data(iris, package="Rdimtools")
#' X   = as.matrix(iris[,1:4])
#' lab = as.factor(as.vector(iris[,5]))
#' ids = sample(1:150, 30)
#'
#' Xold = X[setdiff(1:150,ids),]  # 80% of data for training
#' Xnew = X[ids,]                 # 20% of data for testing
#'
#' ## run PCA for train data & use the info for prediction
#' training = do.pca(Xold,ndim=2)
#' Yold     = training$Y
#' Ynew     = Xnew%*%training$projection
#' Yplab    = lab[ids]
#'
#' ## perform out-of-sample prediction
#' Yoos  = oos.linproj(Xold, Yold, Xnew)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(Ynew, pch=19, col=Yplab, main="true prediction")
#' plot(Yoos, pch=19, col=Yplab, main="OOS prediction")
#' par(opar)
#' }
#'
#' @author Kisung You
#' @rdname oos_LINPROJ
#' @export
oos.linproj <- function(Xold, Yold, Xnew){
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

  # # 2. trfinfo
  # if (!is.list(trfinfo)){
  #   stop(paste("* oos.",methodname," : 'trfinfo' should be provided a list.",sep=""))
  # }
  # if ((!('algtype'%in%names(trfinfo)))||(!('type'%in%names(trfinfo)))||(!('mean'%in%names(trfinfo)))||(!('multiplier'%in%names(trfinfo)))){
  #   stop(paste("* oos.",methodname," : 'trfinfo' is an invalid one. Use 'info' output you acquired from one of functions in the package."))
  # }
  # #------------------------------------------------------------------------
  # # COMPUTE : PREPROCESSING
  # # old ones
  # X = aux.oospreprocess(Xold, trfinfo)
  # Y = Yold
  # # new ones
  # XX = aux.oospreprocess(Xnew, trfinfo)
  #
  # #------------------------------------------------------------------------
  # # COMPUTE : MAIN PART FOR LINEAR PROJECTION
  # # 1. X\beta = Y
  # betahat = aux.pinv(X)%*%Y
  # # 2. use the projection matrix for this
  # YY      = XX%*%betahat
  #
  #
  # #------------------------------------------------------------------------
  # ## RETURN !
  # result = list()
  # result$Ynew = YY;
  # return(result)

  #------------------------------------------------------------------------
  ## COMPUTE
  Ynew = oos_linproj(Xold, Yold, Xnew)
  return(Ynew)
}


