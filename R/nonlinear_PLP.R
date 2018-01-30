#' Piecewise Laplacian-based Projection (PLP)
#'
#' \code{do.plp} is an implementation of Piecewise Laplacian-based Projection (PLP) that
#' adopts two-stage reduction scheme with local approximation.
#'
#' First step is to select \eqn{\sqrt{n}} number of control points using \eqn{k}-means algorithm.
#' After selecting control points that play similar roles as representatives of the entire data points,
#' it performs classical multidimensional scaling.
#'
#' For the rest of the data other than control points,
#' Laplacian Eigenmaps (\code{\link{do.lapeig}}) is then applied to high-dimensional data points
#' lying in neighborhoods of each control point. Embedded low-dimensional local manifold is then
#' aligned to match their coordinates as of their counterparts from classical MDS.
#'
#' @section Notes:
#' \emph{Random Control Points} : The performance of embedding using PLP heavily relies on
#' selection of control points, which is contingent on the performance of \eqn{k}-means
#' clustering.
#'
#' \emph{User Interruption} : PLP is actually an interactive algorithm that a user should be able to intervene
#' intermittently. Such functionality is, however, sacrificed in this version.
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null" and three options of "center", "decorrelate", or "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#'
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate default dataset
#' X <- aux.gensamples()
#'
#' ## run 5 times under the same setting
#' par(mfrow=c(1,5))
#' for (i in 1:5){
#'     out = do.plp(X,ndim=2,type=c("proportion",0.2))
#'     pm  = paste("iteration ",i,sep="")
#'     plot(out$Y[,1],out$Y[,2],main=pm)
#' }
#' }
#'
#' @references
#' \insertRef{paulovich_piece_2011}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_PLP
#' @export
do.plp <- function(X,ndim=2,preprocess="null",type=c("proportion",0.20)){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  if ((!is.numeric(ndim))||(ndim<2)||(ndim>=ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.plp : 'ndim' is a positive integer in [2,#(covariates)).")
  }

  k = as.integer(ndim)
  n = nrow(X)      # number of data points
  d = ncol(X)      # original data dimension
  m = round(sqrt(n))  # number of subset points

  # 2. Parameters
  # 2-1. Common
  #   preprocess     : 'null'(default),'center','whiten','decorrelate'
  # 2-2. plp-LE
  #   type        : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #                 set default for "proportion" of 20%

  if (!is.element(preprocess,c("null","center","whiten","decorrelate"))){
    stop("* do.plp : 'preprocess' should have one of 4 values.")
  }
  nbdtype = type

  # 3. Run
  #   3-1. preprocess
  if (preprocess=="null"){
    trfinfo = list()
    trfinfo$type = "null"
    pX = as.matrix(X,nrow=nrow(X))
  } else {
    tmplist = aux.preprocess(X,type=preprocess)
    trfinfo = tmplist$info
    pX      = tmplist$pX
  }
  #   3-2. sample selection : no random sampling : I will use kmeans
  kmeans  = kmeans(pX,m)
  cluster = kmeans$cluster # assignment information

  #   3-3. control points
  CP = array(0,c(1,d))
  CPvec = c() # class information
  CPidx = c() # corresponding indices
  for (i in 1:m){
    idxi = which(cluster==i)
    samplevec = sample(idxi,round(sqrt(length(idxi))))
    CP = rbind(CP,pX[samplevec,])
    CPvec = c(CPvec,rep(i,times=round(sqrt(length(idxi)))))
    CPidx = c(CPidx,samplevec)
  }
  CP = CP[-1,]

  #   3-4. dimension reduction : MDS
  outputMDS = do.mds(CP,ndim=k,preprocess="center")
  Ytmp  = aux.preprocess(outputMDS$Y,type="center")
  Ycontrol  = Ytmp$pX

  #   3-5. Piecewise Laplacian Projection
  Youtput = array(0,c(n,k))
  #   3-5-extra : pseudoinverse
  pseudoinv = function(X, tol = sqrt(.Machine$double.eps))
  {
    ## Generalized Inverse of a Matrix
    dnx <- dimnames(X)
    if(is.null(dnx)) dnx <- vector("list", 2)
    s <- svd(X)
    nz <- s$d > tol * s$d[1]
    structure(
      if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
      dimnames = dnx[2:1])
  }
  for (i in 1:m){
    # 3-5-1. indexing
    idxall = which(cluster==i)
    idxContrl  = CPidx[which(CPvec==i)]
    idxOthers   = setdiff(idxall,idxContrl)

    # 3-5-2. partial data
    pdata = rbind(pX[idxContrl,],pX[idxOthers,]) # <- this needs to be LEigen into k dims
    pctrl = Ycontrol[which(CPvec==i),]

    # 3-5-3. apply Laplacian Eigenmaps
    outputLE = do.lapeig(pdata,ndim=k,type=nbdtype,symmetric="union")
    pdataY   = outputLE$Y

    # 3-5-4. align
    meanA = colMeans(pdataY[1:length(idxContrl),])
    meanB = colMeans(pctrl)

    cA = pdataY - matrix(rep(meanA,times=nrow(pdata)),nrow=nrow(pdata),byrow=TRUE)
    cB = pctrl - matrix(rep(meanB,times=nrow(pctrl)),nrow=nrow(pctrl),byrow=TRUE)
    C = pseudoinv(cA[1:length(idxContrl),]) %*% cB

    # 3-5-5. assign
    Youtput[idxall,] = (cA %*% C) + matrix(rep(meanB,times=nrow(cA)),nrow=nrow(cA),byrow=TRUE)
  }

  # 4. return output
  result = list()
  result$Y = Youtput
  trfinfo$algtype = "nonlinear"
  result$trfinfo  = trfinfo
  return(result)
}
