#' Independent Component Analysis
#'
#' \code{do.ica} is an R implementation of FastICA algorithm, which aims at
#' finding weight vectors that maximize a measure of non-Gaussianity of projected data.
#' FastICA is initiated with pre-whitening of the data. Single and multiple component
#' extraction are both supported. For more detailed information on ICA and FastICA algorithm,
#' see this \href{https://en.wikipedia.org/wiki/FastICA}{Wikipedia} page.
#'
#' In most of ICA literature, we have \deqn{S = X*W} where \eqn{W} is an unmixing matrix for
#' the given data \eqn{X}. In order to preserve consistency throughout our package, we changed
#' the notation; \code{Y} a projected matrix for \eqn{S} and \code{projection} for unmixing matrix \eqn{W}.
#'
#' @param X an \code{(n-by-p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type nonquadratic function, one of \code{"logcosh"},\code{"exp"}, or \code{"poly"} be chosen.
#' @param tpar a numeric parameter for \code{logcosh} and \code{exp} parameters that should be close to 1.
#' @param sym a logical value; \code{FALSE} for not using symmetric decorrelation, \code{TRUE} otherwise.
#' @param tol stopping criterion for iterative update.
#' @param redundancy a logical value; \code{TRUE} for removing \code{NA} values after prewhitening, \code{FALSE} otherwise.
#' @param maxiter maximum number of iterations allowed.
#'
#' \describe{
#' \item{Y}{an \code{(n-by-ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \code{(p-by-ndim)} whose columns are basis for projection.}
#' }
#'
#'@examples
#'# generate data
#'## in order to pass CRAN pretest, n is set to be small.
#'X <- aux.gensamples(n=28)
#'
#'## 1. use logcosh function for transformation
#'output1 <- do.ica(X,ndim=2,type="logcosh")
#'
#'## 2. use exponential function for transformation
#'output2 <- do.ica(X,ndim=2,type="exp")
#'
#'## 3. use polynomial function for transformation
#'output3 <- do.ica(X,ndim=2,type="poly")
#'
#'## Visualize three different projections
#'par(mfrow=c(1,3))
#'plot(output1$Y[,1],output1$Y[,2],main="logcosh")
#'plot(output2$Y[,1],output2$Y[,2],main="exp")
#'plot(output3$Y[,1],output3$Y[,2],main="poly")
#'
#' @references
#' \insertRef{hyvarinen_independent_2001}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_ICA
#' @export
do.ica <- function(X,ndim=2,type="logcosh",tpar=1,sym=FALSE,tol=1e-6,redundancy=TRUE,maxiter=100){
  # For this, prewhitening is default
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)

  # 2. preprocessing
  tmplist = aux.preprocess(X,type="whiten")
  trfinfo = tmplist$info
  trfinfo$algtype = "linear"
  pX      = tmplist$pX
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||(is.na(ndim))||(is.infinite(ndim))){
    stop("* do.ica : ndim is invalid - should be in [1,#(covariates))")
  }
  ndim = as.integer(min(max(1,ndim),ncol(X)))
  # 3. Parameter Setup
  #   3-1. type : logcosh(1), exp(2), poly(3)
  if (type=="logcosh"){tnum = as.integer(1)}
  else if (type=="exp"){tnum = as.integer(2)}
  else if (type=="poly"){tnum = as.integer(3)}
  else {stop("* do.ica : 'type' parameter should be either 'logcosh','exp', or 'poly'.")}

  #   3-2. tpar : parameter for logcosh or exp
  if (!is.numeric(tpar)||is.na(tpar)||is.infinite(tpar)){
    stop("* do.ica : 'tpar' should be a positive real number around 1.")
  }
  if (tnum==1){
    if ((tpar<1)||(tpar>2)){
      stop("* do.ica : for 'logcosh' type, 'tpar' should be in [1,2]")
    }
  } else if (tnum==2){
    if ((tpar<=0)||(tpar>=2)){
      stop("* do.ica : for 'exp' type, 'tpar' should be (0,2).")
    }
  }

  #   3-3. sym : symmetric decorrelation (default is FALSE)
  if (!is.logical(sym)){
    stop("* do.ica : 'sym' should be a logical parameter")
  }

  #   3-4.tol : tolerance level in iteration
  mepsil = .Machine$double.eps
  if (is.na(tol)||is.infinite(tol)||(tol<mepsil)||(tol>=1)||(!is.numeric(tol))){
    stop("* do.ica : 'tol' should be in [machine epsilon,1)")
  }
  tol = as.double(min(mepsil*1e+6,tol))
  #   3-5. redundancy : remove NA after prewhitening (TRUE)
  if (!is.logical(redundancy)){
    stop("* do.ica : 'redundancy' is a logical indicator.")
  }
  #   3-6. maxiter : 100
  if ((!is.numeric(maxiter))||(maxiter<5)||(is.na(maxiter))||(is.infinite(maxiter))){
    stop("* do.ica : 'maxiter' is a positive integer greater than or equal to 5.")
  }
  maxiter = as.integer(maxiter)


  # 4. Main Code : check NA values
  tpX = t(pX)
  if (redundancy){
    rmNAtpX = tpX[rowSums(is.na(tpX))==0,]
    ndim = min(ndim,nrow(rmNAtpX))
    output = method_ica(rmNAtpX,ndim,maxiter,tol,tnum,tpar,sym)
  } else {
    output = method_ica(tpX,ndim,maxiter,tol,tnum,tpar,sym)
  }

  # 5. Result : S(M*ndim) = X(M*p)*W(p*ndim)
  #             Y         = X*projection
  result = list()
  result$Y = t(output$S)
  result$trfinfo = trfinfo
  result$projection = (output$W)
  return(result)
}
