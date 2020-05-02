#' Nearest Neighbor Projection
#'
#' Nearest Neighbor Projection is an iterative method for visualizing high-dimensional dataset
#' in that a data is sequentially located in the low-dimensional space by maintaining
#' the triangular distance spread of target data with its two nearest neighbors in the high-dimensional space.
#' We extended the original method to be applied for arbitrarily low-dimensional space. Due the generalization,
#' we opted for a global optimization method of \emph{Differential Evolution} (\code{\link[RcppDE]{DEoptim}}) within in that it may add computational burden to certain degrees.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' ## let's compare with other methods
#' out1 <- do.nnp(X, ndim=2)      # NNP
#' out2 <- do.pca(X, ndim=2)      # PCA
#' out3 <- do.lamp(X, ndim=2)     # LAMP
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=label, main="NNP")
#' plot(out2$Y, pch=19, col=label, main="PCA")
#' plot(out3$Y, pch=19, col=label, main="LAMP")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{tejada_improved_2003}{Rdimtools}
#'
#' @rdname nonlinear_NNP
#' @author Kisung You
#' @concept nonlinear_methods
#' @export
do.nnp <- function(X, ndim=2, preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  ########################################################################
  ## 1. Type Checking
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.nnp : 'ndim' is a positive integer in [1,#(covariates)).")}
  algpreprocess = match.arg(preprocess)

  ########################################################################
  ## 2. Preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  D       = as.matrix(dist(pX))

  ########################################################################
  ## 3. Initial Projection
  included = rep(FALSE, n)
  pY       = array(0,c(n,ndim))
  initid   = aux.findmaxidx(D)
  initQ    = as.integer(initid[1])  # two farthest points
  initR    = as.integer(initid[2])

  included[initid] = TRUE # adjust included ones
  pY[initQ,] = c(D[initQ,initR]/2, rep(0,(ndim-1)))
  pY[initR,] = c(-D[initQ,initR]/2, rep(0,(ndim-1)))

  toberun = setdiff(1:n, initid)


  # for general case, I'm using this.
  minfunc <- function(z,center1,center2,rad1,rad2){
    return(abs(sqrt(sum((z-center1)^2))-rad1)+abs(sqrt(sum((z-center2)^2))-rad2))
  }
  Dctrl   = RcppDE::DEoptim.control(trace=FALSE)
  bdlower = rep(-D[initQ,initR], ndim)
  bdupper = rep(D[initQ,initR], ndim)
  ########################################################################
  ## 4. Main Iteration
  for (i in 1:length(toberun)){
    # 4-1. current index
    idx = toberun[i]
    # 4-2. find two nearest ones
    currentincluded = which(included) # index for currently available ones
    if (length(currentincluded)==2){
      idq = currentincluded[1]
      idr = currentincluded[2]
    } else {
      currentdist     = D[idx,currentincluded]
      bottom2         = currentdist[order(currentdist)[1:3]]
      idq = currentincluded[(currentdist==bottom2[2])]
      idr = currentincluded[(currentdist==bottom2[3])]
      if (length(idq)>1){
        idq = as.integer(idq[1])
      }
      if (length(idr)>1){
        idr = as.integer(idr[1])
      }
    }

    # 4-3. get ready for distances
    dxq = as.double(D[idx,idq])
    dxr = as.double(D[idx,idr])
    dxsum = (dxq+dxr)
    dqrlow = sqrt(sum(as.vector(pY[idq,]-pY[idr,])^2))
    # 4-4. type branching
    q = as.vector(pY[idq,])
    r = as.vector(pY[idr,])
    if (dxsum==dqrlow){
      pY[idx,] = ((r-q)*dxq/(dxq+dxr))+q
    } else if (dxsum>dqrlow){
      runDEoptim = RcppDE::DEoptim(minfunc, lower=bdlower, upper=bdupper, control=Dctrl, center1=q, center2=r, rad1=dxq, rad2=dxr)
      pY[idx,] = as.vector((runDEoptim$optim$bestmem))
    } else { # now, we have two cases
      if (dxq<dxr){
        pY[idx,] = q+((r-q)*0.5)
      } else {
        pY[idx,] = r+((q-r)*0.5)
      }
    }
    # 4-5. update index
    included[idx] = TRUE
  }

  ########################################################################
  ## 5. return output
  result = list()
  result$Y = pY
  trfinfo$algtype = "nonlinear"
  result$trfinfo  = trfinfo
  return(result)
}

## Force Algorithm is not applied.

