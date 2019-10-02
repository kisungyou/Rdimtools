#' Interactive Document Map
#'
#' Interactive Document Map originates from text analysis to generate maps of documents by placing
#' similar documents in the same neighborhood. After defining pairwise distance with cosine similarity,
#' authors asserted to use either \code{NNP} or \code{FastMap} as an engine behind.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param engine either \code{NNP} or \code{FastMap}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @seealso \code{\link{do.nnp}}, \code{\link{do.fastmap}}
#'
#' @examples
#' \dontrun{
#' ## load iris data
#' data(iris)
#' X <- as.matrix(iris[,1:4])
#'
#' ## let's compare with other methods
#' out1 <- do.pca(X, ndim=2)
#' out2 <- do.sne(X, ndim=2)
#' out3 <- do.idmap(X, ndim=2, engine="NNP")
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="PCA")
#' plot(out2$Y[,1], out2$Y[,2], main="SNE")
#' plot(out3$Y[,1], out3$Y[,2], main="IDMAP")
#' }
#'
#' @references
#' \insertRef{minghim_content-based_2006}{Rdimtools}
#' @export
do.idmap <- function(X, ndim=2, preprocess=c("null","center","scale","cscale","whiten","decorrelate"), engine=c("NNP","FastMap")){
  ########################################################################
  ## 1. Type Checking
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.idmap : 'ndim' is a positive integer in [1,#(covariates)).")}
  algpreprocess = match.arg(preprocess)

  ########################################################################
  ## 2. Preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  ########################################################################
  ## 3. Preliminary Computation of Cosine Similarity
  D = array(0,c(n,n))
  for (i in 1:n){
    D[i,i] = 0
  }
  for (i in 1:(n-1)){
    ti = as.vector(pX[i,])
    for (j in (i+1):n){
      tj = as.vector(pX[j,])

      scos   = sum(ti*tj)/sqrt(sum(ti^2)*sum(tj^2))
      theval = sqrt(2*(1-scos))
      D[i,j] = theval
      D[j,i] = theval
    }
  }

  ########################################################################
  ## 4. typechecking
  engine = match.arg(engine)
  if (engine=="FastMap"){
    k=ndim
    Dold=D
    Dnew    = array(0,c(n,n))
    output  = array(0,c(n,ndim))
    for (i in 1:k){
      # 3-1. find row and column index for maximal element
      maxidx = aux.findmaxidx(Dold)
      ida    = maxidx[1]
      idb    = maxidx[2]
      # 3-2. precompute some values
      dab2  = (sum(as.vector(pX[ida,]-pX[idb,])^2))
      if (dab2 > (sqrt(123*.Machine$double.eps))){
        dab   = sqrt(dab2)
        # 3-3. compute coefficient
        for (j in 1:n){
          output[j,i] = (sum(as.vector(pX[ida,]-pX[j,])^2) + dab2 - sum(as.vector(pX[idb,]-pX[j,])^2))/(2*dab)
        }
      } # or, leave it as zero

      # 3-4. update D : compute and alter
      for (it1 in 1:n){
        for (it2 in 1:n){
          Dnew[it1,it2] = sqrt((Dold[it1,it2]^2) - ((output[it1,i]-output[it2,i])^2))
        }
      }
      Dold = Dnew
    }
    result = list()
    result$Y = output
    trfinfo$algtype = "nonlinear"
    result$trfinfo  = trfinfo
  } else if (engine=="NNP"){ #######-----------------------################
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
    result = list()
    result$Y = pY
    trfinfo$algtype = "nonlinear"
    result$trfinfo  = trfinfo
  }

  ########################################################################
  ## 5. return output
  return(result)
}
