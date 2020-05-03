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
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(nn1$mask); title("knn with k=20")
#' image(nn2$mask); title("enn with radius=1")
#' image(nn3$mask); title("proportion of ratio=0.4")
#' par(opar)
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
    if (all(method=="minkowski")){
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
  if (all(type[1]=="knn")){
    flag = TRUE
    nnpar = as.double(type[2])
    if (!is.numeric(nnpar)||is.na(nnpar)||is.infinite(nnpar)||(nnpar<1)||(nnpar>ndata)){
      stop("* aux.graphnbd : for knn method, k should be a number in [1,#(data)]")
    } else {
      nnpar = as.integer(nnpar)
    }
  } else if (all(type[1]=="enn")){
    flag = FALSE
    nnpar = as.double(type[2])
    if (!is.numeric(nnpar)||is.na(nnpar)||is.infinite(nnpar)||(nnpar<=0)){
      stop("* aux.graphnbd : for enn method, radius should be a positive real number.")
    } else {
      nnpar = as.double(nnpar)
    }
  } else if (all(type[1]=="proportion")){
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
  if (all(symmetric=="union")){
    outMask = matrix(as.logical(outMask + t(outMask)),nrow=ndata)
  } else if (all(symmetric=="intersect")){
    outMask = matrix(as.logical(outMask * t(outMask)),nrow=ndata)
  } else if (all(symmetric=="asymmetric")){
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
  } else if (all(type[1]=="enn")){
    flag = FALSE
    nnpar = as.double(type[2])
    if (!is.numeric(nnpar)||is.na(nnpar)||is.infinite(nnpar)||(nnpar<=0)){
      stop("* aux.graphnbd : for enn method, radius should be a positive real number.")
    } else {
      nnpar = as.double(nnpar)
    }
  } else if (all(type[1]=="proportion")){
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
  if (all(symmetric=="union")){
    outMask = matrix(as.logical(outMask + t(outMask)),nrow=ndata)
  } else if (all(symmetric=="intersect")){
    outMask = matrix(as.logical(outMask * t(outMask)),nrow=ndata)
  } else if (all(symmetric=="asymmetric")){
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

