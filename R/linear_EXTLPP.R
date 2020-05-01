#' Extended Locality Preserving Projection
#'
#' Extended Locality Preserving Projection (EXTLPP) is an unsupervised
#' dimension reduction algorithm with a bit of flavor in adopting
#' discriminative idea by nature. It raises a question on the data points
#' at \emph{moderate} distance in that a Z-shaped function is introduced in
#' defining similarity derived from Euclidean distance.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param numk the number of neighboring points for k-nn graph construction.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## generate data
#' set.seed(100)
#' X <- aux.gensamples(n=75)
#'
#' ## run Extended LPP with different neighborhood graph
#' out1 <- do.extlpp(X, numk=5)
#' out2 <- do.extlpp(X, numk=10)
#' out3 <- do.extlpp(X, numk=25)
#'
#' ## Visualize three different projections
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="EXTLPP::k=5")
#' plot(out2$Y, main="EXTLPP::k=10")
#' plot(out3$Y, main="EXTLPP::k=25")
#' par(opar)
#'
#' @references
#' \insertRef{shikkenawis_improving_2012}{Rdimtools}
#'
#' @seealso \code{\link{do.lpp}}
#' @author Kisung You
#' @rdname linear_EXTLPP
#' @concept linear_methods
#' @export
do.extlpp <- function(X, ndim=2, numk=max(ceiling(nrow(X)/10),2),
                      preprocess=c("center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.extlpp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. numk
  numk = as.integer(numk)
  if (!check_NumMM(numk,1,n/2,compact=FALSE)){stop("* do.extlpp : 'numk' should be an integer in [2,nrow(X)/2).")}
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. K-Means Clustering
  kclust     = stats::kmeans(pX, numk)
  clustlabel = kclust$cluster
  clustidx   = list() # for each label, find the corresponding # length-'numk' list
  for (i in 1:numk){
    clustidx[[i]] = which(clustlabel==unique(clustlabel)[i])
  }

  #   3. pairwise distance
  PD = as.matrix(dist(pX))
  vecb = rep(0,numk)
  for (i in 1:numk){
    tgtidx = clustidx[[i]]
    vecb[i] = max(PD[tgtidx,tgtidx])
  }
  veca = rep(min(vecb)/20,numk)

  #   4. compute S
  S = array(0,c(n,n))
  for (i in 1:numk){
    tgtidx = clustidx[[i]]
    tgtmat = as.matrix(PD[tgtidx,tgtidx],nrow=length(tgtidx))
    S[tgtidx,tgtidx] = method_trfextlpp(tgtmat,as.double(veca[i]),as.double(vecb[i]))
  }
  diag(S) = 0.0

  #   5. graph laplaciana and generalized eigenvalue problem
  D = diag(rowSums(S))
  L = D-S

  LHS = t(pX)%*%L%*%pX
  RHS = t(pX)%*%D%*%pX

  #   6. compute Projection Matrix : use lowest ones
  projection = aux.geigen(LHS, RHS, ndim, maximal=FALSE)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}

