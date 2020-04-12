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
# 18. aux.oosprocess       : data processign for oos prediction
# 19. aux.findmaxidx       : find the row and column index of maximal elements
# 20. aux.randpartition    : given 1:n, divide it into K random partitions without replacement
# 21. aux.which.mink       : returns index of smallest
#     aux.which.maxk
# 22. aux.traceratio       : solve trace ratio problem with 2012 Ngo's algorithm
# 23. aux.2scatter         : compute LDA type within- and between- scatter matrices
# 24. aux.subsetid         : generate 'k' subset id

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
#' Preprocessing the data
#'
#' \code{aux.preprocess} can perform one of following operations; \code{"center"}, \code{"scale"},
#' \code{"cscale"}, \code{"decorrelate"} and \code{"whiten"}. See below for more details.
#'
#' @section
#' Operations: we have following operations,
#'  \describe{
#'  \item{\code{"center"}}{subtracts mean of each column so that every variable has mean \eqn{0}.}
#'  \item{\code{"scale"}}{turns each column corresponding to variable have variance \eqn{1}.}
#'  \item{\code{"cscale"}}{combines \code{"center"} and \code{"scale"}.}
#'  \item{\code{"decorrelate"}}{\code{"center"} and sets its covariance term having diagonal entries only.}
#'  \item{\code{"whiten"}}{\code{"decorrelate"} and sets all diagonal elements be \eqn{1}.}
#'  }
#'
#' @param data an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param type one of \code{"center"}, \code{"scale"}, \code{"cscale"}, \code{"decorrelate"} or \code{"whiten"}.
#' @return named list containing:
#' \describe{
#' \item{pX}{an \eqn{(n\times p)} matrix after preprocessing in accordance with \code{type} parameter}
#' \item{info}{a list containing \itemize{
#' \item \code{type:} name of preprocessing procedure.
#' \item \code{mean:} a mean vector of length \eqn{p}.
#' \item \code{multiplier:} a \eqn{(p\times p)} matrix or 1 for "center".}}
#' }
#'
#' @examples
#' \dontrun{
#' ## Generate data
#' X = aux.gensamples()
#'
#' ## 5 types of preprocessing
#' X_center = aux.preprocess(X)
#' X_scale  = aux.preprocess(X,type="scale")
#' X_cscale = aux.preprocess(X,type="cscale")
#' X_decorr = aux.preprocess(X,type="decorrelate")
#' X_whiten = aux.preprocess(X,type="whiten")
#'
#' ## Check with Covariance matrix
#' par(mfrow=c(2,3))
#' image(cov(X)[,3:1],          zlim=c(-5,5)); title("original covariance")
#' image(cov(X_center$pX)[,3:1],zlim=c(-5,5)); title("opt::center")
#' image(cov(X_scale$pX)[,3:1], zlim=c(-5,5)); title("opt::scale")
#' image(cov(X_cscale$pX)[,3:1],zlim=c(-5,5)); title("opt::cscale")
#' image(cov(X_decorr$pX)[,3:1],zlim=c(-5,5)); title("opt::decorrelate")
#' image(cov(X_whiten$pX)[,3:1],zlim=c(-5,5)); title("opt::whiten")
#' }
#'
#' @rdname aux_preprocess
#' @author Kisung You
#' @export
aux.preprocess <- function(data,type=c("center","scale","cscale","decorrelate","whiten")){
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
  type = match.arg(type)

  ## two added methods : 'scale' and 'cscale'
  if (type=="scale"){ # add 1 : "scale"
    tmpdata = t(matinput)
    p       = ncol(tmpdata)

    multiplier = rep(0,p)
    for (i in 1:p){
      multiplier[i] = 1/stats::sd(as.vector(tmpdata[,i]))
    }
    info = list()
    info$type = "scale"
    info$mean = rep(0,p)
    info$multiplier = diag(multiplier)

    result = list()
    result$pX = tmpdata%*%diag(multiplier)
    result$info = info
    return(result)
  } else if (type=="cscale"){ # add 2 : "center" + "scale"
    tmpdata = t(matinput)
    p       = ncol(tmpdata)

    infomean = as.double(colMeans(tmpdata))
    infomultiplier = rep(0,p)
    for (i in 1:p){
      infomultiplier[i] = 1/stats::sd(as.vector(tmpdata[,i]))
    }
    info = list()
    info$type = "cscale"
    info$mean = infomean
    info$multiplier = diag(infomultiplier)

    outdata = array(0,c(nrow(tmpdata),p)) # mean substraction
    for (i in 1:nrow(tmpdata)){
      outdata[i,] = as.vector(tmpdata[i,])-infomean
    }

    result = list()
    result$pX   = (outdata%*%diag(infomultiplier))
    result$info = info
    return(result)
  } else {
    # original methods
    # 'center','decorrelate', or 'whiten'
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
}
#' @keywords internal
#' @noRd
aux.preprocess.hidden <- function(data,type=c("null","center","scale","cscale","decorrelate","whiten"),algtype=c("linear","nonlinear")){
  ## minimal preprocessing
  pptype    = match.arg(type)
  ppalgtype = match.arg(algtype)

  ## null is mine !
  n = nrow(data)
  p = ncol(data)
  if (type=="null"){
    info = list()
    info$type = "null"
    info$mean = rep(0,p)
    info$multiplier = 1

    result      = list()
    result$pX   = data
    result$info = info
    result$info$algtype = ppalgtype
    return(result)
  } else {
    result = aux.preprocess(data,type=pptype)
    result$info$algtype = ppalgtype
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
#'
#' @param dname name of a predefined shape. Should be one of \describe{
#' \item{\code{"swiss"}}{swiss roll}
#' \item{\code{"crown"}}{crown}
#' \item{\code{"helix"}}{helix}
#' \item{\code{"saddle"}}{manifold near saddle point}
#' \item{\code{"ribbon"}}{ribbon}
#' \item{\code{"bswiss"}}{broken swiss}
#' \item{\code{"cswiss"}}{cut swiss}
#' \item{\code{"twinpeaks"}}{two peaks}
#' \item{\code{"sinusoid"}}{sinusoid on the circle}
#' \item{\code{"mobius"}}{mobius strip embedded in \eqn{\mathbf{R}^3}}
#' \item{\code{"R12in72"}}{12-dimensional manifold in \eqn{\mathbf{R}^{12}}}
#' }
#' @param ... extra parameters for the followings #' \tabular{lll}{
#' parameter \tab dname \tab description \cr
#' \code{ntwist} \tab \code{"mobius"} \tab number of twists
#' }
#'
#' @return an \eqn{(n\times 3)} matrix of generated data by row.
#'
#' @references
#' \insertRef{hein_intrinsic_2005}{Rdimtools}
#'
#' \insertRef{vandermaaten_learning_2009}{Rdimtools}
#'
#' @author Kisung You
#' @rdname aux_gensamples
#' @export
aux.gensamples <- function(n=496,noise=0.01,
                           dname=c("swiss","crown","helix","saddle","ribbon","bswiss","cswiss","twinpeaks","sinusoid",
                                   "mobius","R12in72"), ...){
  #params = as.list(environment())
  n     = as.integer(n)
  noise = as.double(noise)
  dtype = (match.arg(dname))

  ## extra parameters
  extrapar = list(...)
  #   1. ntwist
  if ("ntwist" %in% names(extrapar)){
    ntwist = as.integer(extrapar$ntwist)
  } else {
    ntwist= as.integer(10)
  }
  if ((ntwist<1)||(is.na(ntwist))||(length(ntwist)>1)||(is.infinite(ntwist))){
    stop("* aux.gensamples : 'ntwist' should be a positive integer.")
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
  } else if (dtype=="sinusoid"){
    # hein audibert
    tt = runif(n, min=0, max=(2*3.141592))

    hx = sin(tt)
    hy = cos(tt)
    hz = sin(150*tt)/10

    highD = cbind(hx,hy,hz)+matrix(rnorm(n*3,sd=noise),nrow=n)
  } else if (dtype=="mobius"){
    u = runif(n,-1,1)
    v = runif(n,0,2*3.141592)

    hx = (1+(u/2)*cos(ntwist*v/2))*(cos(v))
    hy = (1+(u/2)*cos(ntwist*v/2))*(sin(v))
    hz = (u/2)*sin(ntwist*v/2)

    highD = cbind(hx,hy,hz)+matrix(rnorm(n*3,sd=noise),nrow=n)
  } else if (dtype=="R12in72"){
    xi = array(0,c(n,24)) # without noise
    for (it in 1:n){
      alpha = runif(12)
      for (i in 1:11){
        xi[it,(2*i-1)] = alpha[i+1]*cos(2*pi*alpha[i])
        xi[it,2*i]     = alpha[i+1]*sin(2*pi*alpha[i])
      }
      xi[it,23] = alpha[1]*cos(2*pi*alpha[12])
      xi[it,24] = alpha[1]*sin(2*pi*alpha[12])
    }

    highD = cbind(xi,xi,xi) + matrix(rnorm(n*72,sd=noise),nrow=n)
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
  npoints = as.integer(npoints)
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
  dnormval = base::norm(diag(p)-Pid, type = "F")
  if (dnormval > 1e-10){
    return(qr.Q(qr(P)))
  } else {
    return(P)
  }
}
#' @keywords internal
#' @noRd
aux.adjprojection <- function(P){
  PP = aux.adjqr(P)
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
  # new 1. use 'maotai's mimic function : increasing order as well
  gfun  = getFromNamespace("hidden_geigen","maotai")
  geigs = gfun(top, bottom, normalize=TRUE)
  nobj  = length(geigs$values)
  ndim  = round(ndim)

  # new 2. separate eigenvectors accordingly to the orders
  if (maximal){
    vectors = geigs$vectors[,nobj:(nobj-ndim+1)]
  } else {
    vectors = geigs$vectors[,1:ndim]
  }

  # new 3. return
  if (ndim==1){
    return(matrix(vectors, ncol=1))
  } else {
    return(vectors)
  }

  # # 1. first, run CPP with RcppArmadillo -> change to 'geigen' package
  # # geigs = aux_geigen(top, bottom) # Armadillo goes Decreasing order
  # geigs = geigen::geigen(top, bottom) # geigen goes Increasing order
  # maxp  = length(geigs$values)
  #
  # # 2. separate values and vectors; change correspondingly for Decreasing order
  # values  = geigs$values[maxp:1]
  # vectors = geigs$vectors[,maxp:1]
  #
  # if (maximal==TRUE){
  #   partvals = values[1:ndim]
  # } else {
  #   partvals = values[maxp:(maxp-ndim+1)]
  # }
  # if (all(base::abs(base::Im(partvals))<(100*(.Machine$double.eps)))){
  #   values  = base::Re(values)
  #   vectors = aux.adjprojection(base::Re(vectors))
  # } else {
  #   stop("* aux.geigen : generalized eigenvalue problem returned imaginary eigenvalues.")
  # }
  #
  #
  # # 3. branching case
  # if (ndim > 1){
  #   if (maximal==TRUE){
  #     projection = vectors[,1:ndim]
  #   } else {
  #     projection = vectors[,maxp:(maxp-ndim+1)]
  #   }
  # } else {
  #   if (maximal==TRUE){
  #     vecsol = vectors[,1]
  #     projection = matrix(vecsol/sqrt(sum(vecsol*vecsol)))
  #   } else {
  #     vecsol = vectors[,maxp]
  #     projection = matrix(vecsol/sqrt(sum(vecsol*vecsol)))
  #   }
  # }
  # return(projection)
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
  m = round(ndim)
  r = round(aux_rank(B)) # as.integer(Matrix::rankMatrix (B))

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

  A = matrix(A,nrow=nrow(A))
  if (is.vector(B)){
    B = matrix(B,ncol=1)
  } else {
    B = matrix(B,nrow=nrow(B))
  }
  preconditioner = matrix(preconditioner,nrow=nrow(preconditioner))
  sparseflag = FALSE

  # sparseformats = c("dgCMatrix","dtCMatrix","dsCMatrix")
  # if ((class(A)%in%sparseformats)||(class(B)%in%sparseformats)||(class(preconditioner)%in%sparseformats)){
  #   A = Matrix(A,sparse=TRUE)
  #   B = Matrix(B,sparse=TRUE)
  #   preconditioner = Matrix(preconditioner,sparse=TRUE)
  #   sparseflag = TRUE
  # } else {
  #   A = matrix(A,nrow=nrow(A))
  #   if (is.vector(B)){
  #     B = matrix(B)
  #   } else {
  #     B = matrix(B,nrow=nrow(B))
  #   }
  #   preconditioner = matrix(preconditioner,nrow=nrow(preconditioner))
  #   sparseflag = FALSE
  # }
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
      vecB = as.vector(B[,i])
      tmpres = linsolve.bicgstab.single(A,vecB,xinit,reltol,maxiter,preconditioner)
      # if (!sparseflag){
      #
      # } else {
      #   vecB = Matrix (B[,i],sparse=TRUE)
      #   tmpres = linsolve.bicgstab.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner)
      # }
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

# 18. aux.oospreprocess ---------------------------------------------------
#     data processing for out-of-sample prediction
#' @keywords internal
#' @noRd
aux.oospreprocess <- function(data, trfinfo){
  ## 0. parameter
  n = nrow(data)
  p = ncol(data)
  output = array(0,c(n,p))
  ## 1. extract mean
  meanvec = as.vector(trfinfo$mean)
  for (i in 1:n){
    output[i,] = (as.vector(data[i,])-meanvec)
  }
  ## 2. multiplier
  multiplier = trfinfo$multiplier
  if (is.matrix(multiplier)){
    output = output%*%multiplier
  } else {
    output = output*multiplier
  }
  return(output)
}

# 19. aux_findmaxidx       : find the row and column index of maxi --------
#' @keywords internal
#' @noRd
aux.findmaxidx <- function(A){
  # 19. aux_findmaxidx       : find the row and column index of maximal elements
  output = which(A==max(A), arr.ind=TRUE)
  if (nrow(output)>1){
    return(output[1,])
  } else {
    return(output)
  }
}


# 20. aux.randpartition ---------------------------------------------------
#     given 1:n, divide it into K random partitions without replacement
#' @keywords internal
#' @noRd
aux.randpartition <- function(n, K){
  output = list()
  if (K==1){
    output = list()
    output[[1]] = 1:n
  } else {
    listall = 1:n
    singleK = round(n/K)
    for (i in 1:(K-1)){
      output[[i]] = sample(listall, singleK, replace = FALSE)
      listall = setdiff(listall, output[[i]])
    }
    output[[K]] = listall
  }
  return(output)
}


# 21. aux.which.mink ------------------------------------------------------
#     aux.which.maxk
#' @keywords internal
#' @noRd
aux.which.mink <- function(x, k=1){
  return(order(x)[1:k])
}
#' @keywords internal
#' @noRd
aux.which.maxk <- function(x, k=1){
  return(order(x,decreasing = TRUE)[1:k])
}


# 22. aux.traceratio  : solve trace ratio problem with 2012 Ngo's  --------
#' @keywords internal
#' @noRd
aux.traceratio <- function(A, B, dim, eps, maxiter){
  ## in the the language
  n = nrow(A)
  p = dim

  ## prepare the initializer
  Vold = qr.Q(qr(matrix(rnorm(n*p),ncol=p)))
  rhoold = 0
  for (i in 1:maxiter){
    Vnew   = RSpectra::eigs(A-rhoold*B,p,which="LR")$vectors
    rhonew = sum(diag(t(Vnew)%*%A%*%Vnew))/sum(diag(t(Vnew)%*%B%*%Vnew))

    rhoinc = abs(rhonew-rhoold)
    Vold   = Vnew
    rhoold = rhonew

    if (rhoinc < eps){
      break
    }
  }

  ## let's try to return !
  return(Vold)
}


# 23. aux.2scatter --------------------------------------------------------
#     data should be provided as a matrix (columns are variables)
#' @keywords internal
#' @noRd
aux.2scatter <- function(pX, label){
  # 0. extra information
  if (is.vector(pX)){
    n = length(pX)
    p = 1
  } else {
    n = nrow(pX)
    p = ncol(pX)
  }


  # 1. extract label information
  label   = round(label)
  ulabel  = unique(label)
  datlist = list()
  for (i in 1:length(ulabel)){
    if (is.vector(pX)){
      datlist[[i]] = pX[(label==ulabel[i])]
    } else {
      datlist[[i]] = pX[(label==ulabel[i]),]
    }
  }
  # 2. compute two types of scatter matrices
  #   2-1. error matrix/E/error variance
  scattermat <- function(x){
    return(cov(x)*(nrow(x)-1))
  }
  matE = array(0,c(p,p))
  for (i in 1:length(ulabel)){
    if (is.vector(datlist[[i]])){
      tgt = as.matrix(datlist[[i]])
    } else {
      tgt = datlist[[i]]
    }
    matE = matE + cov(tgt)*(nrow(tgt)-1)
  }
  matH = array(0,c(p,p))
  if (is.vector(datlist[[1]])){
    meanlist = lapply(datlist, base::mean)
  } else {
    meanlist = lapply(datlist, colMeans)
  }
  if (is.vector(pX)){
    meantott = base::mean(pX)
  } else {
    meantott = colMeans(pX)
  }
  for (i in 1:length(ulabel)){
    meandiff = as.vector(meanlist[[i]]-meantott)
    if (is.vector(datlist[[i]])){
      matH = matH + length(datlist[[i]])*outer(meandiff,meandiff)
    } else {
      matH = matH + nrow(datlist[[i]])*outer(meandiff,meandiff)
    }
  }

  output = list()
  if (is.vector(pX)){
    output$within = as.double(matE)
    output$between = as.double(matH)
  } else {
    output$within = matE
    output$between = matH
  }

  return(output)
}


# 24. aux.subsetid --------------------------------------------------------
#' @keywords internal
aux.subsetid <- function(n, k){
  x = sample(1:n)
  return(split(x, sort(x%%k)))
}
