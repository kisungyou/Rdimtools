#' Diffusion Maps
#'
#' \code{do.dm} discovers low-dimensional manifold structure embedded in high-dimensional
#' data space using Diffusion Maps (DM). It exploits diffusion process and distances in data space to find
#' equivalent representations in low-dimensional space.
#'
#' @param X an \code{(n-by-p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is ``null'' and three options of ``center'',``decorrelate'', or ``whiten''
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param diffscale a scaling parameter for diffusion kernel. Default is 1 and should be a positive real number with 0 inclusive.
#' @param timescale a target scale whose value represents behavior of heat kernels at time \emph{t}. Default is 1 and should be a positive real number >0.
#' @param threshold a numerical threshold for making laplacian sparse. Default is 1e-6.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \code{(n-by-ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{eigvals}{a vector of eigenvalues for Markov transition matrix.}
#' }
#'
#'
#'@examples
#'## generate Swiss Roll data of 28 data points.
#'## in order to pass CRAN pretest, n is set to be small.
#'X <- aux.gensamples(n=28)
#'
#'## Various combinations of (diffscale, timescale)
#'out1 <- do.dm(X,ndim=2,diffscale=1,timescale=0.1)
#'out2 <- do.dm(X,ndim=2,diffscale=5,timescale=0.1)
#'out3 <- do.dm(X,ndim=2,diffscale=10,timescale=0.1)
#'out4 <- do.dm(X,ndim=2,diffscale=1,timescale=1)
#'out5 <- do.dm(X,ndim=2,diffscale=5,timescale=1)
#'out6 <- do.dm(X,ndim=2,diffscale=10,timescale=1)
#'out7 <- do.dm(X,ndim=2,diffscale=1,timescale=10)
#'out8 <- do.dm(X,ndim=2,diffscale=5,timescale=10)
#'out9 <- do.dm(X,ndim=2,diffscale=10,timescale=10)
#'
#'## Visualize 9 different combinations
#'par(mfrow=c(3,3))
#'plot(out1$Y[,1],out1$Y[,2],main="(diff=1,time=0.1)")
#'plot(out2$Y[,1],out2$Y[,2],main="(diff=5,time=0.1)")
#'plot(out3$Y[,1],out3$Y[,2],main="(diff=10,time=0.1)")
#'plot(out4$Y[,1],out4$Y[,2],main="(diff=1,time=1)")
#'plot(out5$Y[,1],out5$Y[,2],main="(diff=5,time=1)")
#'plot(out6$Y[,1],out6$Y[,2],main="(diff=10,time=1)")
#'plot(out7$Y[,1],out7$Y[,2],main="(diff=1,time=10)")
#'plot(out8$Y[,1],out8$Y[,2],main="(diff=5,time=10)")
#'plot(out9$Y[,1],out9$Y[,2],main="(diff=10,time=10)")
#'
#'@references
#'\insertRef{nadler_diffusion_2005}{Rdimtools}
#'
#'\insertRef{coifman_diffusion_2006}{Rdimtools}
#'
#'
#' @rdname nonlinear_DM
#' @author Kisung You
#' @export
do.dm <- function(X,ndim=2,preprocess="null",diffscale=1,timescale=1,threshold=1e-7){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.dm : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  k = as.integer(ndim)
  n = nrow(X)
  d = ncol(X)

  # 2. Parameters
  # 2-1. Common
  #   preprocess     : 'null'(default),'center','whiten','decorrelate'
  # 2-2. Diffusion Maps only
  #   diffscale      : 1(default) or a real number >= 0
  #   timescale      : 1(default) or a real number >  0
  #   threshold      : 1e-7(default)

  if (!is.element(preprocess,c("null","whiten","center","decorrelate"))){
    stop("* do.dm : 'preprocess' should have one of 4 values.")
  }
  if (preprocess=="null"){
    preprocess = FALSE
  }
  if (!is.numeric(diffscale)|(diffscale<0)|is.infinite(diffscale)){
    stop("* do.dm : 'diffscale' should be a real number >= 0.")
  }
  if (!is.numeric(timescale)|(timescale<=0)|is.infinite(timescale)){
    stop("* do.dm : 'timescale' should be a positive real number > 0.")
  }
  if (!is.numeric(threshold)|is.infinite(threshold)|(threshold<0)){
    stop("* do.dm : 'threshold' value should be a small number >= 0.")
  }
  threshold = max((min(1e-7,threshold)),(.Machine$double.eps)*10)


  # 3. Run
  #   3-1. preprocess
  if (preprocess==FALSE){
    trfinfo = list()
    trfinfo$type = "null"
    pX = as.matrix(X,nrow=nrow(X))
  } else {
    tmplist = aux.preprocess(X,type=preprocess)
    trfinfo = tmplist$info
    pX      = tmplist$pX
  }
  #   3-2. diffusion matrix
  K = exp(-((as.matrix(dist(pX)))^2)/diffscale)
  #   3-3. I'll use symmetric version
  Dalpha = diag(1/sqrt(rowSums(K)))
  A = Dalpha %*% K %*% Dalpha
  A[(abs(A)<threshold)] = 0
  #   3-4. SVD decomposition
  rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
  }
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  svdA = svd(A,nu=(k+1),nv=(k+1))
  svdmultiplier = rep.col(svdA$u[,1],(k+1))
  psi  = (svdA$u)/svdmultiplier # right eigenvectors of Markov matrix
  phi  = (svdA$u)*svdmultiplier # left  eigenvectors of Markov matrix
  eigvals = svdA$d
  #   3-5. compute embedding
  lambdat = ((svdA$d[2:(k+1)])^timescale)
  lambdat = rep.row(lambdat,nrow(psi))
  Y       = (psi[,-1] * lambdat)

  # 4. output
  result = list()
  result$Y = Y
  trfinfo$algtype   = "nonlinear"
  result$trfinfo = trfinfo
  result$eigvals = eigvals
  return(result)
}
