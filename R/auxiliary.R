# LIST OF AUXILIARY FUNCTIONS
# 00. aux.typecheck        : check whether the data is poorly given
# 01. aux.preprocess       : preprocess of centering, decorrelation or whitening
# 02. aux.gensamples       : generate a few popular data examples
# 03. aux.graphnbd         : construct nearest-neighborhood graph
#     aux.graphnbdD        : distance matrix is already given
# 04. aux.shortestpath     : compute shortest paths given neighborhood graph
# 05. aux.MaxMinLandmark   : choose a single landmark point
# 06. aux.kernelcov        : build K and centerd K matrix for kernel tricks
# 07. aux.eigendec         : use Armadillo in a descending order
# 08. aux.pkgstat          : show number of available functions.
# 09. aux.kernelcentering  : centering the kernel/gram matrix
# 10. aux.kernelprojection : given uncentered gram matrix, find the projected data
#                            note that it results (ndim-by-N) matrix, columns are projected vectors.
# 11. aux.adjprojection    : adjust projection matrix by simply normalizing each column
#     aux.adjqr            : qr decomposition is used for adjusting
# 12. aux.nbdlogical       : find homogeneous and heterogeneous neighborhood indexing
# 13. aux.geigen           : geigen in my taste
# 14. aux.featureindicator : generate (p-by-ndim) indicator matrix for projection
# 15. aux.traceratio.max   : compute trace ratio problem for maximal basis
# 16. aux.pinv             : use SVD and NumPy scheme
# 17. aux.bicgstab         : due to my stupidity, now Rlinsolve can't be used.

#  ------------------------------------------------------------------------
# 0. AUX.TYPECHECK
#  ------------------------------------------------------------------------
#' @noRd
#' @keywords internal
aux.typecheck <- function(data, verbose=FALSE){
  # data frame into matrix
  matinput = as.matrix(data)
  if ((dim(matinput)[1]==1)||(dim(matinput)[2]==1)){
    if (verbose==TRUE){
      message("WARNING : input data should be matrix, not a vector.")
    }
    matinput = matrix(matinput)
  }

  # check if there exists any NA, Inf, -Inf
  if (!is.numeric(matinput)){
    warning("ERROR : input data should be numeric arrays.")
    return(FALSE)
  }
  if (any(is.infinite(matinput))||any(is.na(matinput))){
    warning("ERROR : input data should contain no Inf or NA values.")
    return(FALSE)
  }
  return(TRUE)
}


#  ------------------------------------------------------------------------
# 1. AUX.PREPROCESS
#  ------------------------------------------------------------------------
#' Centering, decorrelating, or whitening of the data
#'
#' \code{aux.preprocess} can perform one of three popular operations; centering, decorrelating,
#'  and whitening of data. See below for more details.
#'
#' \code{"center"} option subtracts mean of each column so that every variable has
#' mean 0. After centering, option \code{"decorrelate"} sets the data matrix
#' to have diagonal covariance terms only. \code{"whiten"} option sets the sample
#' covariance to have all diagonal terms equal to 1, equally weighting each variable.
#'
#' @param data an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param type one of "center", "decorrelate", or "whiten". See below for more details.
#' @return named list containing:
#' \describe{
#' \item{pX}{an \eqn{(n\times p)} matrix after preprocessing in accordance with \code{type} parameter}
#' \item{info}{a list containing \itemize{
#' \item \code{type:} name of preprocessing procedure.
#' \item \code{mean:} a mean vector of length \eqn{p}.
#' \item \code{multiplier:} a \eqn{(p\times p)} matrix for "decorrelate" or "whiten" or 1 for "center".}}
#' }
#'
#' @examples
#' \dontrun{
#' ## Generate data
#' X = aux.gensamples()
#'
#' ## 3 types of preprocessing
#' X_center = aux.preprocess(X)
#' X_decorr = aux.preprocess(X,type="decorrelate")
#' X_whiten = aux.preprocess(X,type="whiten")
#'
#' ## Check with Covariance matrix
#' par(mfrow=c(1,3))
#' image(cov(X_center$pX),zlim=c(-50,50)); title("center");
#' image(cov(X_decorr$pX),zlim=c(-50,50)); title("decorrelate");
#' image(cov(X_whiten$pX),zlim=c(-50,50)); title("whitening")
#' }
#'
#' @rdname aux_preprocess
#' @author Kisung You
#' @export
aux.preprocess <- function(data,type="center"){
  # data : (n-by-d)
  #   ARMA with CPP
  #   input & output are (d-by-n)
  #   Don't forget to transpose at the end.
  if (is.data.frame(data)){
    matinput = t(as.matrix(data))
  } else if (is.matrix(data)){
    matinput = t(data)
  } else {
    warning("WARNING : input should be either dataframe or matrix.")
  }

  # 'center','decorrelate', or 'whiten'
  if (!is.element(type,c("center","decorrelate","whiten"))){
    stop("* aux.preprocess: 'type' should be one of 3 options; center, decorrelate, or whiten.")
  }
  matoutput = tryCatch(
    {
      if (type=="center"){
        aux_preprocess(matinput,as.integer(1))
      } else if (type=="decorrelate"){
        aux_preprocess(matinput,as.integer(2))
      } else {
        aux_preprocess(matinput,as.integer(3))
      }
    }, error=function(cond){
      return(NA)
    }, warning=function(cond){
      return(NA)
    }
  )
  if (length(matoutput)==1){
    if (is.na(matoutput)){
      result = NA
      return(result)
    }
  } else {
    # now we have 2 lists
    info   = list()
    info$type = matoutput$type
    info$mean = matoutput$mean
    info$multiplier = matoutput$multiplier

    result = list()
    result$pX = t(matoutput$output)
    result$info = info
    return(result)
  }
}


#  ------------------------------------------------------------------------
# 2. AUX.GENSAMPLES
#  ------------------------------------------------------------------------
#' Generate model-based samples
#'
#' It generates samples from predefined 3-d shapes, set by \code{dname} parameter.
#' Also incorporated a functionality to add white noise with degree \code{noise}.
#'
#' @param n the number of points to be generated.
#' @param noise level of additive white noise.
#' @param dname name of a predefined shape. Should be one of \describe{
#' \item{\code{"swiss"}}{swiss roll}
#' \item{\code{"crown"}}{crown}
#' \item{\code{"helix"}}{helix}
#' \item{\code{"saddle"}}{manifold near saddle point}
#' \item{\code{"ribbon"}}{ribbon}
#' \item{\code{"bswiss"}}{broken swiss}
#' \item{\code{"cswiss"}}{cut swiss}
#' \item{\code{"twinpeaks"}}{two peaks}
#' }
#' @return an \eqn{(n\times 3)} matrix of generated data by row.
#' @examples
#' ## Generate samples for three different shapes
#' d1 = aux.gensamples(dname="twinpeaks",noise=0.01)
#' d2 = aux.gensamples(dname="ribbon",noise=0.01)
#' d3 = aux.gensamples(dname="crown", noise=0.01)
#'
#' \dontrun{
#' casenames = c("swiss","crown","helix","saddle","ribbon","bswiss","cswiss","twinpeaks")
#' for (i in 1:length(casenames)){
#'   data = aux.gensamples(n=sample(1000:2000,1),noise=runif(1)[1],dname=casenames[i])
#' }
#' }
#' @references
#' \insertRef{van_der_maaten_dimensionality_2009}{Rdimtools}
#'
#'
#' @author Kisung You
#' @rdname aux_gensamples
#' @export
aux.gensamples <- function(n=496,noise=0.1,dname="swiss"){
  params = as.list(environment())
  dtype = params$dname
  casenames = c("swiss","crown","helix","saddle","ribbon","bswiss","cswiss","twinpeaks")
  if (!is.element(dtype,casenames)){
    message("function gensamples:: use one of specificed types")
    stop()
  }

  # 2. type branching
  if (dtype=="swiss"){
    lowx = ((3*pi)/2)*(1+2*runif(n))

    highx = lowx*cos(lowx)
    highy = lowx*sin(lowx)
    highz = 30*runif(n)

    highD = cbind(highx,highy,highz) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (dtype=="crown"){
    lowD = (2*pi)*(1:n)/n

    highx = (2+sin(8*lowD))*sin(lowD)
    highy = (2+sin(8*lowD))*cos(lowD)
    highz = sin(8*lowD)

    highD = cbind(highx,highy,highz) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (dtype=="helix"){
    lowD = (2*pi)*(1:n)/n

    highx = (2+cos(8*lowD))*cos(lowD)
    highy = (2+cos(8*lowD))*sin(lowD)
    highz = sin(8*lowD)

    highD = cbind(highx,highy,highz) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (dtype=="saddle"){
    lowD = 1-2*matrix(rnorm(2*n),c(n,2))

    hx = lowD[,1]
    hy = lowD[,2]
    highz = sin(pi*hx)*tanh(3*hy)

    highD = cbind(lowD,highz) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (dtype=="ribbon"){
    lowD = (2*pi)*(1:n)/n

    hx = cos(lowD)
    hy = sin(lowD)
    h = 5*runif(n)

    highD = cbind(hx,hx*hy,h) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (dtype=="bswiss"){
    nceil = ceiling(n/2)
    lowD1 = (3*pi/2)*(1+0.8*runif(nceil))
    lowD2 = (2*pi/2)*(1+0.8*runif(n-nceil)+0.6)
    lowD  = c(lowD1,lowD2)

    height = 30*runif(n)

    hx = lowD*cos(lowD)
    hy = height
    hz = lowD*sin(lowD)

    highD = cbind(hx,hy,hz) + matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (dtype=="cswiss"){
    r = array(0,c(1,n))
    for (i in 1:n){
      pass = FALSE
      while (!pass){
        rr = runif(1)
        if (runif(1) > rr){
          r[i] = rr
          pass = TRUE
        }
      }
    }

    lowD = t((3*pi/2)*(1+(2*r)))
    height = 21*runif(n)

    hx = lowD*cos(lowD)
    hy = height
    hz = lowD*sin(lowD)

    highD = cbind(hx,hy,hz)+matrix(rnorm(n*3,sd=noise),c(n,3))
  } else if (dtype=="twinpeaks"){
    xi = runif(n)
    yi = runif(n)

    hx = 1-2*xi
    hy = sin(pi-2*pi*(xi))
    hz = tanh(3-6*yi)

    highD = cbind(hx,hy,hz)+matrix(rnorm(n*3,sd=noise),c(n,3))
  }
  return(highD)
}

#  ------------------------------------------------------------------------
# 3. AUX.GRAPHNBD
#  ------------------------------------------------------------------------
#' Construct Nearest-Neighborhood Graph
#'
#' Given data, it first computes pairwise distance (\code{method}) using one of measures
#' defined from \code{\link[stats]{dist}} function. Then, \code{type} controls how nearest neighborhood
#' graph should be constructed. Finally, \code{symmetric} parameter controls how
#' nearest neighborhood graph should be symmetrized.
#'
#' @section Nearest Neighbor(NN) search:
#' Our package supports three ways of defining nearest neighborhood. First is
#' \emph{knn}, which finds \code{k} nearest points and flag them as neighbors.
#' Second is \emph{enn} - epsilon nearest neighbor - that connects all the
#' data poinst within a certain radius. Finally, \emph{proportion} flag is to
#' connect proportion-amount of data points sequentially from the nearest to farthest.
#'
#' @section Symmetrization:
#' In many graph setting, it starts from dealing with undirected graphs.
#' NN search, however, does not necessarily guarantee if symmetric connectivity
#' would appear or not. There are two easy options for symmetrization;
#' \code{intersect} for connecting two nodes if both of them are
#' nearest neighbors of each other and \code{union} for only either of them to be present.
#'
#'
#' @param data an \eqn{(n\times p)} data matrix.
#' @param method type of distance to be used. See also \code{\link[stats]{dist}}.
#' @param type a defining pattern of neighborhood criterion. One of \describe{
#' \item{c("knn", k)}{knn with \code{k} a positive integer.}
#' \item{c("enn", radius)}{enn with a positive radius.}
#' \item{c("proportion", ratio)}{takes an \code{ratio} in (0,1) portion of edges to be connected.}
#' }
#' @param symmetric either ``intersect'' or ``union'' for symmetrization, or ``asymmetric''.
#' @param pval a \eqn{p}-norm option for Minkowski distance.
#' @return a named list containing \describe{
#' \item{mask}{a binary matrix of indicating existence of an edge for each element.}
#' \item{dist}{corresponding distance matrix. \code{-Inf} is returned for non-connecting edges.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Generate data
#' X = aux.gensamples()
#'
#' ## Test three different types of neighborhood connectivity
#' nn1 = aux.graphnbd(X,type=c("knn",20))         # knn with k=20
#' nn2 = aux.graphnbd(X,type=c("enn",1))          # enn with radius = 1
#' nn3 = aux.graphnbd(X,type=c("proportion",0.4)) # connecting 40% of edges
#'
#' ## Visualize
#' par(mfrow=c(1,3))
#' image(nn1$mask); title("knn with k=20")
#' image(nn2$mask); title("enn with radius=1")
#' image(nn3$mask); title("proportion of ratio=0.4")
#' }
#'
#' @author Kisung You
#' @rdname aux_graphnbd
#' @export
aux.graphnbd <- function(data,method="euclidean",type=c("proportion",0.1),symmetric="union",pval=2.0){
  # 1. type check & ...
  aux.typecheck(data)
  if ((!is.numeric(pval))||is.na(pval)||is.infinite(pval)||(pval<=0)){
    stop("* aux.graphnbd : pval should be a positive number as norm sign.")
  }

  # 2. compute distance matrix D
  if (is.element(method,c("euclidean","maximum","manhattan","canberra","binary","minkowski"))){
    if (method=="minkowski"){
      D = as.matrix(dist(data,method="minkowski",p=pval))
    } else {
      D = as.matrix(dist(data,method))
    }
  } else {
    stop("* aux.graphnbd : 'distance' parameter is not valid. Please refer to ?dist")
  }

  # 3.type : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #   note that in case of ratio, it can also be transformed into k
  #   so that we need only two-case branching for the later use.
  #   use a flag  : TRUE for knn style / FALSE for enn style
  #         nnpar : k or epsilon
  ndata = nrow(D)
  if (type[1]=="knn"){
    flag = TRUE
    nnpar = as.double(type[2])
    if (!is.numeric(nnpar)||is.na(nnpar)||is.infinite(nnpar)||(nnpar<1)||(nnpar>ndata)){
      stop("* aux.graphnbd : for knn method, k should be a number in [1,#(data)]")
    } else {
      nnpar = as.integer(nnpar)
    }
  } else if (type[1]=="enn"){
    flag = FALSE
    nnpar = as.double(type[2])
    if (!is.numeric(nnpar)||is.na(nnpar)||is.infinite(nnpar)||(nnpar<=0)){
      stop("* aux.graphnbd : for enn method, radius should be a positive real number.")
    } else {
      nnpar = as.double(nnpar)
    }
  } else if (type[1]=="proportion"){
    flag = TRUE
    nnpar = max(as.double(type[2])*ndata,1)
    if (!is.numeric(nnpar)||is.na(nnpar)||is.infinite(nnpar)||(nnpar<1)||(nnpar>ndata)){
      stop("* aux.graphnbd : for 'proportion' method, ratiovalue should be a number in (0,1]")
    } else {
      nnpar = as.integer(min((ceiling(nnpar)),ndata))
    }
  } else {
    stop("* aux.graphnbd : type pair is not well defined. This must be in a correct form.")
  }

  outMask = array(0,c(ndata,ndata))
  if (flag){
    for (i in 1:ndata){
      tgt = D[i,]
      idx = as.vector(which(tgt<=max((sort(tgt,decreasing=FALSE))[1:min((nnpar+1),ndata)]),arr.ind=T))
      outMask[i,idx] = 1
    }
  } else {
    for (i in 1:ndata){
      tgt = D[i,]
      idx = which(tgt < nnpar)
      outMask[i,idx] = 1
    }
  }
  for (i in 1:ndata){
    outMask[i,i] = 0
  }
  outMask = matrix(as.logical(as.integer(outMask)),nrow=ndata)

  # 4. symmetric="union"'intersect','union', or 'asymmetric'
  if (symmetric=="union"){
    outMask = matrix(as.logical(outMask + t(outMask)),nrow=ndata)
  } else if (symmetric=="intersect"){
    outMask = matrix(as.logical(outMask * t(outMask)),nrow=ndata)
  } else if (symmetric=="asymmetric"){
    outMask = outMask
  } else {
    stop("* aux.graphnbd : 'symmetric' option is invalid.")
  }

  if (sum(!outMask)==0){
    message("* aux.graphnbd : this graph has no connecting edges.")
  }

  # 5. Computation : -Inf means they are not included
  outD1 = outMask * D
  outD2 = (!outMask) * array(-Inf,c(ndata,ndata))
  outD2[which(is.na(outD2))] = 0

  outD = outD1 + outD2
  for (i in 1:ndata){
    outD[i,i] = 0
  }

  result = list()
  result$mask = outMask
  result$dist = outD
  return(result)
}


#' @keywords internal
#' @noRd
aux.graphnbdD <- function(D,type=c("proportion",0.1),symmetric="union",pval=2.0){
  # 1. type check & check parameter 'pval'
  aux.typecheck(D)
  if ((!is.numeric(pval))||is.na(pval)||is.infinite(pval)||(pval<=0)){
    stop("* aux.graphnbd : pval should be a positive number as norm sign.")
  }

  # 2. compute distance matrix D : no I can pass it
  # 3.type : vector of c("knn",k), c("enn",radius), or c("proportion",ratio)
  #   note that in case of ratio, it can also be transformed into k
  #   so that we need only two-case branching for the later use.
  #   use a flag  : TRUE for knn style / FALSE for enn style
  #         nnpar : k or epsilon
  ndata = nrow(D)
  if (type[1]=="knn"){
    flag = TRUE
    nnpar = as.double(type[2])
    if (!is.numeric(nnpar)||is.na(nnpar)||is.infinite(nnpar)||(nnpar<1)||(nnpar>ndata)){
      stop("* aux.graphnbd : for knn method, k should be a number in [1,#(data)]")
    } else {
      nnpar = as.integer(nnpar)
    }
  } else if (type[1]=="enn"){
    flag = FALSE
    nnpar = as.double(type[2])
    if (!is.numeric(nnpar)||is.na(nnpar)||is.infinite(nnpar)||(nnpar<=0)){
      stop("* aux.graphnbd : for enn method, radius should be a positive real number.")
    } else {
      nnpar = as.double(nnpar)
    }
  } else if (type[1]=="proportion"){
    flag = TRUE
    nnpar = max(as.double(type[2])*ndata,1)
    if (!is.numeric(nnpar)||is.na(nnpar)||is.infinite(nnpar)||(nnpar<1)||(nnpar>ndata)){
      stop("* aux.graphnbd : for 'proportion' method, ratiovalue should be a number in (0,1]")
    } else {
      nnpar = as.integer(min((ceiling(nnpar)),ndata))
    }
  } else {
    stop("* aux.graphnbd : type pair is not well defined. This must be in a correct form.")
  }

  outMask = array(0,c(ndata,ndata))
  if (flag){
    for (i in 1:ndata){
      tgt = D[i,]
      idx = as.vector(which(tgt<=max((sort(tgt,decreasing=FALSE))[1:min((nnpar+1),ndata)]),arr.ind=T))
      outMask[i,idx] = 1
    }
  } else {
    for (i in 1:ndata){
      tgt = D[i,]
      idx = which(tgt < nnpar)
      outMask[i,idx] = 1
    }
  }
  for (i in 1:ndata){
    outMask[i,i] = 0
  }
  outMask = matrix(as.logical(as.integer(outMask)),nrow=ndata)

  # 4. symmetric="union"'intersect','union', or 'asymmetric'
  if (symmetric=="union"){
    outMask = matrix(as.logical(outMask + t(outMask)),nrow=ndata)
  } else if (symmetric=="intersect"){
    outMask = matrix(as.logical(outMask * t(outMask)),nrow=ndata)
  } else if (symmetric=="asymmetric"){
    outMask = outMask
  } else {
    stop("* aux.graphnbd : 'symmetric' option is invalid.")
  }

  if (sum(!outMask)==0){
    message("* aux.graphnbd : this graph has no connecting edges.")
  }

  # 5. Computation : -Inf means they are not included
  outD1 = outMask * D
  outD2 = (!outMask) * array(-Inf,c(ndata,ndata))
  outD2[which(is.na(outD2))] = 0

  outD = outD1 + outD2
  for (i in 1:ndata){
    outD[i,i] = 0
  }

  result = list()
  result$mask = outMask
  result$dist = outD
  return(result)
}

#  ------------------------------------------------------------------------
# 4. AUX.SHORTESTPATH
#  ------------------------------------------------------------------------
# either (n*n) matrix or 'dist' object
#' Find shortest path using Floyd-Warshall algorithm
#'
#' This is a fast implementation of Floyd-Warshall algorithm to find the
#' shortest path in a pairwise sense using 'RcppArmadillo'. A logical input
#' is also accepted.
#'
#' @param dist either an \eqn{(n\times n)} matrix or a \code{dist} class object.
#' @return an \eqn{(n\times n)} matrix containing pairwise shortest path.
#' @examples
#' \dontrun{
#' ## Generate 10-sample data
#' X = aux.gensamples(n=10)
#'
#' ## Find knn graph with k=3
#' Xgraph = aux.graphnbd(X,type=c("knn",3))
#'
#' ## Separately use binarized and real distance matrices
#' W1 = aux.shortestpath(Xgraph$mask)
#' W2 = aux.shortestpath(Xgraph$dist)
#'
#' par(mfrow=c(1,2))
#' image(W1); title("from binarized")
#' image(W2); title("from Euclidean distance")
#' }
#'
#' @author Kisung You
#' @references Floyd, R.W. (1962) \emph{Algorithm 97: Shortest Path}. Commincations of the ACMS, Vol.5(6):345.
#' @rdname aux_shortestpath
#' @export
aux.shortestpath <- function(dist){
  # class determination
  if (class(dist)=="matrix"){
    distnaive = dist
  } else if (class(dist)=="dist"){
    distnaive = as.matrix(dist)
  } else {
    stop("* aux.shortestpath : input 'dist' should be either (n*n) matrix or 'dist' class object.")
  }
  # consider logical input
  if (any(is.logical(distnaive))){
    distnaive = distnaive*1
  }
  # set as -Inf for 0 values
  mepsil  = .Machine$double.eps
  distnaive[which(distnaive<5*mepsil)] = -Inf
  distgeo   = aux_shortestpath(distnaive)
  return(distgeo)
}


#  ------------------------------------------------------------------------
# 5. AUX.MAXMINLANDMARK
#  ------------------------------------------------------------------------
#' @noRd
#' @keywords internal
aux.MaxMinLandmark <- function(X,npoints,pdflag=FALSE){
  # 5-1. setting
  if (nrow(X)<=npoints){
    stop("ERROR : npoints should be smaller than the total number of original data points.")
  }

  # 5-2. initialize with two starting points
  landmarkidx = array(0,c(1,npoints))
  nX = nrow(X)
  seqnp = seq_len(nX)
  idx12 = sample(1:nX,2)

  landmarkidx[1:2] = idx12     # recorded random points
  seqnp = setdiff(seqnp,idx12) # sequence left over

  # 5-3. main computation
  if (pdflag==TRUE){
    pD = X
  } else {
    pD = as.matrix(dist(X))
  }
  if (npoints>2){
    for (it in 3:npoints){
      # 5-3-1. marginalize the data
      plandmark  = landmarkidx[1:(it-1)]
      currentidx = aux_landmarkMaxMin(pD, plandmark, seqnp)

      # 5-3-2. update
      landmarkidx[it] = currentidx
      seqnp = setdiff(seqnp,currentidx)
    }
  }

  # 5-4. return results
  return(landmarkidx)
}

#  ------------------------------------------------------------------------
# 6. AUX.KERNELCOV
#  ------------------------------------------------------------------------
#' Build a centered kernel matrix K
#'
#' From the celebrated Mercer's Theorem, we know that for a mapping \eqn{\phi}, there exists
#' a kernel function - or, symmetric bilinear form, \eqn{K} such that \deqn{K(x,y) = <\phi(x),\phi(y)>} where \eqn{<,>} is
#' standard inner product. \code{aux.kernelcov} is a collection of 20 such positive definite kernel functions, as
#' well as centering of such kernel since covariance requires a mean to be subtracted and
#' a set of transformed values \eqn{\phi(x_i),i=1,2,\dots,n} are not centered after transformation.
#' Since some kernels require parameters - up to 2, its usage will be listed in arguments section.
#'
#' @details
#' There are 20 kernels supported. Belows are the kernels when given two vectors \eqn{x,y}, \eqn{K(x,y)}
#' \describe{
#' \item{linear}{\eqn{=<x,y>+c}}
#' \item{polynomial}{\eqn{=(<x,y>+c)^d}}
#' \item{gaussian}{\eqn{=exp(-c\|x-y\|^2)}, \eqn{c>0}}
#' \item{laplacian}{\eqn{=exp(-c\|x-y\|)}, \eqn{c>0}}
#' \item{anova}{\eqn{=\sum_k exp(-c(x_k-y_k)^2)^d}, \eqn{c>0,d\ge 1}}
#' \item{sigmoid}{\eqn{=tanh(a<x,y>+b)}}
#' \item{rational quadratic}{\eqn{=1-(\|x-y\|^2)/(\|x-y\|^2+c)}}
#' \item{multiquadric}{\eqn{=\sqrt{\|x-y\|^2 + c^2}}}
#' \item{inverse quadric}{\eqn{=1/(\|x-y\|^2+c^2)}}
#' \item{inverse multiquadric}{\eqn{=1/\sqrt{\|x-y\|^2+c^2}}}
#' \item{circular}{\eqn{=
#' \frac{2}{\pi} arccos(-\frac{\|x-y\|}{c}) - \frac{2}{\pi} \frac{\|x-y\|}{c}\sqrt{1-(\|x-y\|/c)^2}
#' }, \eqn{c>0}}
#' \item{spherical}{\eqn{=
#' 1-1.5\frac{\|x-y\|}{c}+0.5(\|x-y\|/c)^3
#' }, \eqn{c>0}}
#' \item{power/triangular}{\eqn{=-\|x-y\|^d}, \eqn{d\ge 1}}
#' \item{log}{\eqn{=-\log (\|x-y\|^d+1)}}
#' \item{spline}{\eqn{=
#' \prod_i (
#' 1+x_i y_i(1+min(x_i,y_i)) - \frac{x_i + y_i}{2} min(x_i,y_i)^2
#' + \frac{min(x_i,y_i)^3}{3}
#' )
#' }}
#' \item{Cauchy}{\eqn{=\frac{c^2}{c^2+\|x-y\|^2}}}
#' \item{Chi-squared}{\eqn{=\sum_i \frac{2x_i y_i}{x_i+y_i}}}
#' \item{histogram intersection}{\eqn{=\sum_i min(x_i,y_i)}}
#' \item{generalized histogram intersection}{\eqn{=sum_i min(
#' |x_i|^c,|y_i|^d
#' )}}
#' \item{generalized Student-t}{\eqn{=1/(1+\|x-y\|^d)}, \eqn{d\ge 1}}
#' }
#'
#' @param X an \eqn{(n\times p)} data matrix
#' @param ktype a vector containing the type of kernel and parameters involved. Below the usage is
#' consistent with description
#' \describe{
#' \item{linear}{\code{c("linear",c)}}
#' \item{polynomial}{\code{c("polynomial",c,d)}}
#' \item{gaussian}{\code{c("gaussian",c)}}
#' \item{laplacian}{\code{c("laplacian",c)}}
#' \item{anova}{\code{c("anova",c,d)}}
#' \item{sigmoid}{\code{c("sigmoid",a,b)}}
#' \item{rational quadratic}{\code{c("rq",c)}}
#' \item{multiquadric}{\code{c("mq",c)}}
#' \item{inverse quadric}{\code{c("iq",c)}}
#' \item{inverse multiquadric}{\code{c("imq",c)}}
#' \item{circular}{\code{c("circular",c)}}
#' \item{spherical}{\code{c("spherical",c)}}
#' \item{power/triangular}{\code{c("power",d)}}
#' \item{log}{\code{c("log",d)}}
#' \item{spline}{\code{c("spline")}}
#' \item{Cauchy}{\code{c("cauchy",c)}}
#' \item{Chi-squared}{\code{c("chisq")}}
#' \item{histogram intersection}{\code{c("histintx")}}
#' \item{generalized histogram intersection}{\code{c("ghistintx",c,d)}}
#' \item{generalized Student-t}{\code{c("t",d)}}
#' }
#'
#' @return a named list containing \describe{
#' \item{K}{a \eqn{(p\times p)} kernelizd gram matrix.}
#' \item{Kcenter}{a \eqn{(p\times p)} centered version of \code{K}.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate data
#' X = aux.gensamples(n=100)
#'
#' ## extra parameters do not matter for no-parameter kernels
#' A1 = aux.kernelcov(X,c("spline"))
#' A2 = aux.kernelcov(X,c("spline",1,2)) # these numbers will be disregarded.
#' print(paste("* aux.kernelcov : abs.diff.",norm(A1$K-A2$K,"f"),"in norm"))
#' }
#'
#'
#' @references Hofmann, T., Scholkopf, B., and Smola, A.J. (2008) \emph{Kernel methods in
#' machine learning}. arXiv:math/0701907.
#' @author Kisung You
#' @rdname aux_kernelcov
#' @export
aux.kernelcov <- function(X,ktype){
  # 6-1. 20 ktype is supported
  if (length(ktype)==1){
    defaultflag = TRUE
  } else {
    defaultflag = FALSE
  }
  if (length(ktype)>3){
    stop("* aux.kernelcov : maximum 2 parameters are used.")
  }
  if (typeof(ktype[1])!="character"){
    stop("* aux.kernelcov : first argument should be a name of kernel chosen.")
  }
  if ((ktype[1]=="circular")&&(ncol(X)!=2)){
    stop("* aux.kernelcov : circular kernel is for 2-d data only.")
  }
  if ((ktype[1]=="spherical")&&(ncol(X)!=3)){
    stop("* aux.kernelcov : spherical kernel is for 3-d data only.")
  }

  kchrname = ktype[1]
  switch(kchrname,
         # 1. "linear" - linear kernel
         linear={
           knumber = as.integer(1)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
         },
         #  2. "polynomial" - polynomial kernel
         polynomial={
           knumber = as.integer(2)
           if (defaultflag==TRUE){
             par1 = 1.0
             par2 = 3
           } else {
             if (length(ktype)==2){
               stop("* aux.kernelcov : 'polynomial' kernel requires two parameters.")
             }
             par1 = as.numeric(ktype[2])
             par2 = as.numeric(ktype[3])
           }
         },
         #  3. "gaussian" - gaussian kernel
         gaussian={
           knumber = as.integer(3)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
             if ((par1 <= 0)||(is.na(par1))||(is.infinite(par1))){
               stop("* aux.kernelcov : 'gaussian' scale parameter should be a positive value.")
             }
           }
           par2 = 0.0
         },
         #  4. "laplacian" - laplacian kernel
         laplacian = {
           knumber = as.integer(4)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
             if ((par1 <= 0)||(is.na(par1))||(is.infinite(par1))){
               stop("* aux.kernelcov : 'gaussian' scale parameter should be a positive value.")
             }
           }
           par2 = 0.0
         },
         #  5. "anova" - ANOVA kernel
         anova = {
           knumber = as.integer(5)
           if (defaultflag==TRUE){
             par1 = 1.0
             par2 = 1.0
           } else {
             if (length(ktype)==2){
               stop("* aux.kernelcov : 'anova' kernel requires two parameters.")
             }
             par1 = as.numeric(ktype[2])
             par2 = as.numeric(ktype[3])
             if ((par1<=0)||(par2<1)){
               stop("* aux.kernelcov : 'anova' parameters are invalid.")
             }
           }
         },
         #  6. "sigmoid"
         sigmoid = {
           knumber = as.integer(6)
           if (defaultflag==TRUE){
             par1 = 1/nrow(X)
             par2 = 1.0
           } else {
             if (length(ktype)==2){
               stop("* aux.kernelcov : 'sigmoid' kernel requires two parameters.")
             }
             par1 = as.numeric(ktype[2])
             par2 = as.numeric(ktype[3])
           }
         },
         #  7. "rq" - rational quadratic
         rq = {
           knumber = as.integer(7)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
         },
         #  8. "mq" - multiquadric
         mq = {
           knumber = as.integer(8)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
         },
         #  9. "iq" - inverse quadric
         iq = {
           knumber = as.integer(9)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
         },
         # 10. "imq" - inverse multiquadric
         imq = {
           knumber = as.integer(10)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
         },
         # 11. "circular" - circular kernel
         #    use par2 as pi = 3.141592
         circular = {
           knumber = as.integer(11)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           if (par1 <= 0){
             stop("* aux.kernelcov : parameter has to be a positive number.")
           }
           par2 = pi
           if (ncol(X)!=2){
             stop("* aux.kernelcov : 'circular' kernel should be used on 2-dimensional data.")
           }
         },
         #  12. "spherical" - spherical kernel
         #      use par2 as pi = 3.141592
         spherical = {
           knumber = as.integer(12)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           if (par1 <= 0){
             stop("* aux.kernelcov : parameter has to be a positive number.")
           }
           par2 = pi
           if (ncol(X)!=3){
             stop("* aux.kernelcov : 'spherical' kernel should be used on 3-dimensional data.")
           }
         },
         #  13. "power" - power/triangular kernel
         power = {
           knumber = as.integer(13)
           if (defaultflag==TRUE){
             par1 = 2.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
           if (par1 < 1){
             stop("* aux.kernelcov : 'power' parameter should be >= 1.")
           }
         },
         #  14. "log" - log kernel
         log = {
           knumber = as.integer(14)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
           if (par1 < 1){
             stop("* aux.kernelcov : 'log' parameter should be >= 1.")
           }
         },
         #  15. "spline" - No Parameter Kernel
         spline = {
           knumber = as.integer(15)
           par1 = 0.0
           par2 = 0.0
         },
         #  16. "cauchy" - cauchy kernel
         cauchy = {
           knumber = as.integer(16)
           if (defaultflag==TRUE){
             par1 = 1.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
           if (par1 < 0){
             stop("* aux.kernelcov : 'cauchy' parameter should be >= 0.")
           }
         },
         #  17. "chisq" - Chi-Squared kernel
         chisq = {
           knumber = as.integer(17)
           par1 = 0.0
           par2 = 0.0
         },
         #  18. "histintx" - histogram intersection
         histintx = {
           knumber = as.integer(18)
           par1 = 0.0
           par2 = 0.0
           if (any((X<0))){
             stop("* aux.kernelcov : 'histintx' can be used for data with positive values only.")
           }
         },
         #  19. "ghistintx" - generalized histogram intersection
         ghistintx = {
           knumber = as.integer(19)
           if (defaultflag==TRUE){
             par1 = 1.0
             par2 = 1.0
           } else {
             if (length(ktype)==2){
               stop("* aux.kernelcov : 'ghistintx' kernel requires two parameters.")
             }
             par1 = as.numeric(ktype[2])
             par2 = as.numeric(ktype[3])
           }
         },
         #  20. "t" - student-t distribution
         t = {
           knumber = as.integer(20)
           if (defaultflag==TRUE){
             par1 = 2.0
           } else {
             par1 = as.numeric(ktype[2])
           }
           par2 = 0.0
           if (par1<1){
             stop("* aux.kernelcov : for student t distribution, we need a parameter >= 1.")
           }
         },
         stop("* aux.kernelcov : invalid kernel type. choose one of three or leave it blank.")
  )
  tX = t(X)
  result = aux_kernelcov(tX,knumber,par1,par2)
  return(result);
}

# 7. eigendecomposition : Armadillo + Descending order --------------------
#' @noRd
#' @keywords internal
aux.eigendec <- function(X){
  if (nrow(X)!=ncol(X)){
    stop("ERROR : a given matrix is not a square form.")
  }
  output = aux_eigendecomposition(X);

  result = list()
  result$eigval = rev(output$eigval)
  result$eigvec = output$eigvec[,rev(seq_len(length(output$eigval)))]
  return(result)
}


# pkgstat : show the number of available functions. ---------------------
#' Show the number of functions for \pkg{Rdimtools}.
#'
#' This function is mainly used for tracking progress for this package.
#'
#' @examples
#' ## run with following command
#' aux.pkgstat()
#'
#' @rdname aux_pkgstat
#' @export
aux.pkgstat <- function(){
  ndo  = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "do."))))
  nest = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "est."))))
  naux = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "aux."))))
  noos = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "oos."))))

  mdo  = paste("*  cat1{do. } manifold learning techniques           : ",ndo,sep="")
  mest = paste("*  cat2{est.} intrinsic dimension estimation methods : ",nest,sep="")
  moos = paste("*  cat3{oos.} out-of-sample projection methods       : ",noos,sep="")
  maux = paste("*  cat4{aux.} auxiliary functions available          : ",naux,sep="")
  print("* Number of functions available in Rdimtools package")
  print(mdo)
  print(mest)
  print(moos)
  print(maux)
}


# 09. aux.kernelcentering -------------------------------------------------
#     centering the kernel/gram matrix
#' @keywords internal
#' @noRd
aux.kernelcentering <- function(K){
  N = nrow(K)
  if (ncol(K)!=N){
    stop("* aux.kernelcentering : an input K should be a square matrix.")
  }
  onesN  = array(1,c(N,N))/N
  Ktilde = K-(onesN%*%K)-(K%*%onesN)+(onesN%*%K%*%onesN)
  return(Ktilde)
}


# 10. aux.kernelprojection ------------------------------------------------
# given uncentered gram matrix, find the projected data
# note that it returns (ndim-by-N)
#' @keywords internal
#' @noRd
aux.kernelprojection <- function(KK, ndim){
  KKcentered = aux.kernelcentering(KK)
  KKceigen   = eigen(KKcentered)
  Y          = (t(KKceigen$vectors[,1:ndim]) %*% KK)
  return(Y)
}


# 11. aux.adjprojection : adjust projection matrix ------------------------
#' @keywords internal
#' @noRd
aux.adjqr <- function(P){
  p = ncol(P)
  Pid = (t(P)%*%P)
  if (max(abs(diag(p)-Pid))>1e-18){
    output = qr.Q(qr(P))
  } else {
    output = P
  }
  return(output)
}
#' @keywords internal
#' @noRd
aux.adjprojection <- function(P){
  n = nrow(P)
  p = ncol(P)
  output = array(0,c(n,p))
  for (i in 1:p){
    cvec = as.vector(P[,i])
    output[,i] = cvec/sqrt(sum(cvec*cvec))
  }
  return(output)
}
# 12. aux.nbdlogical : find homogeneous and heterogeneous neighbor --------
#' @keywords internal
#' @noRd
aux.nbdlogical <- function(X, label, khomo, khet){
  D = as.matrix(dist(X))
  n = nrow(D)
  # 1. homogeneous logical matrix
  logical_hom = array(FALSE,c(n,n))
  for (i in 1:n){
    # 1-1. index of same class
    idxhom = setdiff(which(label==label[i]), i)
    # 1-2. which is the smallest numk ?
    partD = as.vector(D[i,idxhom])
    # 1-3. partially smallest ones
    partidx = which(      partD <= max(sort(partD)[1:max(min(khomo, length(idxhom)),1)])    )
    # 1-4. adjust idxhom
    idxhomadj = idxhom[partidx]
    logical_hom[i,idxhomadj] = TRUE
  }
  # 2. heterogeneous logical matrix
  logical_het = array(FALSE,c(n,n))
  for (i in 1:n){
    idxhet = which(label!=label[i])
    partD = as.vector(D[i, idxhet])
    partidx = which(partD <= max(sort(partD)[1:max(min(khet, length(idxhet)),1)]))
    idxhetadj = idxhet[partidx]
    logical_het[i,idxhetadj] = TRUE
  }
  output = list()
  output$hom = logical_hom
  output$het = logical_het
  return(output)
}

# 13. aux.geigen : now it uses RcppArmadillo ------------------------------
#' @keywords internal
#' @noRd
aux.geigen <- function(top, bottom, ndim, maximal=TRUE){
  # 1. first, run CPP with RcppArmadillo
  geigs = aux_geigen(top, bottom)
  maxp  = length(geigs$values)

  # 2. separate values and vectors
  values  = geigs$values
  vectors = geigs$vectors

  if (maximal==TRUE){
    partvals = values[1:ndim]
  } else {
    partvals = values[maxp:(maxp-ndim+1)]
  }
  if (all(base::abs(base::Im(partvals))<(100*(.Machine$double.eps)))){
    values  = base::Re(values)
    vectors = aux.adjprojection(base::Re(vectors))
  } else {
    stop("* aux.geigen : generalized eigenvalue problem returned imaginary eigenvalues.")
  }


  # 3. branching case
  if (ndim > 1){
    if (maximal==TRUE){
      projection = vectors[,1:ndim]
    } else {
      projection = vectors[,maxp:(maxp-ndim+1)]
    }
  } else {
    if (maximal==TRUE){
      vecsol = vectors[,1]
      projection = matrix(vecsol/sqrt(sum(vecsol*vecsol)))
    } else {
      vecsol = vectors[,maxp]
      projection = matrix(vecsol/sqrt(sum(vecsol*vecsol)))
    }
  }
  return(projection)
}


# 14. aux.featureindicator : generate indicator matrix --------------------
#     generate (p-by-ndim) indicator matrix for projection
#' @keywords internal
#' @noRd
aux.featureindicator <- function(p,ndim,idxvec){
  if (length(idxvec)!=ndim){
    stop("* aux.featureindicator : selection had some problem.")
  }
  output = array(0,c(p,ndim))
  for (i in 1:ndim){
    selectedcolumn = as.integer(idxvec[i])
    output[selectedcolumn,i] = 1
  }
  return(output)
}

# 15. aux.traceratio.max : compute trace ratio problem for maximal --------
#' @keywords internal
#' @noRd
aux.traceratio.max <- function(A,B,ndim,tol=1e-10){
  # -------------------------------------------------------------
  # PRELIMINARY SETTING
  if (nrow(A)!=ncol(A)){
    stop("* aux.traceratio : A is not square.")
  } else {
    d = nrow(A)
  }
  if (nrow(B)!=ncol(B)){
    stop("* aux.traceratio : B is not square")
  }
  if (nrow(B)!=d){
    stop("* aux.traceratio : B is not matching size of A.")
  }
  m = as.integer(ndim)
  r = as.integer(Matrix::rankMatrix(B))

  # -------------------------------------------------------------
  # MAIN BRANCHING
  if (m <= (d-r)){
    Z = RSpectra::eigs(B,(d-r),which="SR")$vectors
    cost = t(Z)%*%A%*%Z
    Zright = RSpectra::eigs(cost, m)$vectors
    output = Z%*%Zright;
  } else {
    lambda1 = sum(diag(A))/sum(diag(B))
    lambda2 = sum(base::eigen(A, only.values=TRUE)$values[1:m])/sum(base::eigen(B, only.values = TRUE)$values[1:m])
    lambda  = ((lambda1+lambda2)/2)
    while ((lambda2 -lambda1)>tol){
      gamma = sum(base::eigen((A-lambda*B), only.values = TRUE)$values[1:m])
      if (gamma > 0){
        lambda1 = lambda
      } else {
        lambda2 = lambda
      }
      lambda = ((lambda1+lambda2)/2)
    }
    output = RSpectra::eigs((A-lambda*B),ndim)$vectors
  }
  return(output)
}

# 16. aux.pinv : use SVD and NumPy Scheme ---------------------------------
# https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Singular_value_decomposition_(SVD)
#' @keywords internal
#' @noRd
aux.pinv <- function(A){
  svdA      = base::svd(A)
  tolerance = (.Machine$double.eps)*max(c(nrow(A),ncol(A)))*as.double(max(svdA$d))

  idxcut    = which(svdA$d <= tolerance)
  invDvec   = (1/svdA$d)
  invDvec[idxcut] = 0

  output = (svdA$v%*%diag(invDvec)%*%t(svdA$u))
  return(output)
}




# 17. aux.bicgstab --------------------------------------------------------
#' @keywords internal
#' @noRd
aux.bicgstab <- function(A,B,xinit=NA,reltol=1e-5,maxiter=1000,
                            preconditioner=diag(ncol(A)),verbose=TRUE){
  ###########################################################################
  # Step 0. Initialization
  if (verbose){
    message("* aux.bicgstab : Initialiszed.")
  }
  if (any(is.na(A))||any(is.infinite(A))||any(is.na(B))||any(is.infinite(B))){
    stop("* aux.bicgstab : no NA or Inf values allowed.")
  }
  sparseformats = c("dgCMatrix","dtCMatrix","dsCMatrix")
  if ((class(A)%in%sparseformats)||(class(B)%in%sparseformats)||(class(preconditioner)%in%sparseformats)){
    A = Matrix(A,sparse=TRUE)
    B = Matrix(B,sparse=TRUE)
    preconditioner = Matrix(preconditioner,sparse=TRUE)
    sparseflag = TRUE
  } else {
    A = matrix(A,nrow=nrow(A))
    if (is.vector(B)){
      B = matrix(B)
    } else {
      B = matrix(B,nrow=nrow(B))
    }
    preconditioner = matrix(preconditioner,nrow=nrow(preconditioner))
    sparseflag = FALSE
  }
  # xinit
  if (is.na(xinit)){
    xinit = matrix(rnorm(ncol(A)))
  } else {
    if (length(xinit)!=ncol(A)){
      stop("* aux.bicgstab : 'xinit' has invalid size.")
    }
    xinit = matrix(xinit)
  }
  ###########################################################################
  # Step 1. Preprocessing
  # 1-1. Neither NA nor Inf allowed.
  if (any(is.infinite(A))||any(is.na(A))||any(is.infinite(B))||any(is.na(B))){
    stop("* aux.bicgstab : no NA, Inf, -Inf values are allowed.")
  }
  # 1-2. Size Argument
  m = nrow(A)
  if (is.vector(B)){
    mB = length(B)
    if (m!=mB){
      stop("* aux.bicgstab : a vector B should have a length of nrow(A).")
    }
  } else {
    mB = nrow(B)
    if (m!=mB){
      stop("* aux.bicgstab : an input matrix B should have the same number of rows from A.")
    }
  }
  if (is.vector(B)){
    B = as.matrix(B)
  }
  # 1-3. Adjusting Case
  if (m > ncol(A)){        ## Case 1. Overdetermined
    B = t(A)%*%B
    A = t(A)%*%A
  } else if (m < ncol(A)){ ## Case 2. Underdetermined
    stop("* aux.bicgstab : underdetermined case is not supported.")
  }
  # 1-4. Preconditioner : only valid for square case
  if (!all.equal(dim(A),dim(preconditioner))){
    stop("* aux.bicgstab : Preconditioner is a size-matching.")
  }
  if (verbose){message("* aux.bicgstab : preprocessing finished ...")}
  ###########################################################################
  # Step 2. Main Computation
  ncolB = ncol(B)
  if (ncolB==1){
    if (!sparseflag){
      vecB = as.vector(B)
      res = linsolve.bicgstab.single(A,vecB,xinit,reltol,maxiter,preconditioner)
    } else {
      vecB = B
      res = linsolve.bicgstab.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner)
    }
  } else {
    x      = array(0,c(ncol(A),ncolB))
    iter   = array(0,c(1,ncolB))
    errors = list()
    for (i in 1:ncolB){
      if (!sparseflag){
        vecB = as.vector(B[,i])
        tmpres = linsolve.bicgstab.single(A,vecB,xinit,reltol,maxiter,preconditioner)
      } else {
        vecB = Matrix(B[,i],sparse=TRUE)
        tmpres = linsolve.bicgstab.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner)
      }
      x[,i]        = tmpres$x
      iter[i]      = tmpres$iter
      errors[[i]]  = tmpres$errors
      if (verbose){
        message(paste("* aux.bicgstab : B's column.",i,"being processed.."))
      }
    }
    res = list("x"=x,"iter"=iter,"errors"=errors)
  }

  ###########################################################################
  # Step 3. Finalize
  if ("flag" %in% names(res)){
    flagval = res$flag
    if (flagval==0){
      if (verbose){
        message("* aux.bicgstab : convergence well achieved.")
      }
    } else if (flagval==1){
      if (verbose){
        message("* aux.bicgstab : convergence not achieved within maxiter.")
      }
    } else {
      if (verbose){
        message("* aux.bicgstab : breakdown.")
      }
    }
    res$flag = NULL
  }
  if (verbose){
    message("* aux.bicgstab : computations finished.")
  }
  return(res)
}
