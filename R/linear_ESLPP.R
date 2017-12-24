#' Extended Supervised Locality Preserving Projection
#'
#' Extended LPP and Supervised LPP are two variants of the celebrated Locality Preserving Projection (LPP) algorithm for dimension
#' reduction. Their combination, Extended Supervised LPP, is a combination of two algorithmic novelties in one that
#' it reflects discriminant information with realistic distance measure via Z-score function.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param numk the number of neighboring points for k-nn graph construction.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center" and other options of "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{eigval}{a vector of eigenvalues corresponding to basis expansion in an ascending order.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' ## generate data of 2 types with clear difference
#' diff = 15
#' dt1  = aux.gensamples(n=123)-diff;
#' dt2  = aux.gensamples(n=123)+diff;
#'
#' ## merge the data and create a label correspondingly
#' Y      = rbind(dt1,dt2)
#' label  = c(rep(1,123), rep(2,123))
#'
#' ## compare LPP, SLPP and ESLPP
#' outLPP  <- do.lpp(Y)
#' outSLPP <- do.slpp(Y, label)
#' outESLPP <- do.eslpp(Y, label)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(outLPP$Y[,1], outLPP$Y[,2], main="LPP")
#' plot(outSLPP$Y[,1], outSLPP$Y[,2], main="SLPP")
#' plot(outESLPP$Y[,1], outESLPP$Y[,2], main="ESLPP")
#' }
#'
#' @references
#' \insertRef{zheng_gabor_2007}{Rdimtools}
#'
#' \insertRef{shikkenawis_improving_2012}{Rdimtools}
#'
#' @seealso \code{\link{do.lpp}}, \code{\link{do.slpp}}, \code{\link{do.extlpp}}
#' @author Kisung You
#' @rdname linear_ESLPP
#' @export
do.eslpp <- function(X, label, ndim=2, numk=max(ceiling(nrow(X)/10),2),
                      preprocess=c("center","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label : check and return a de-factored vector
  #   For this example, there should be no degenerate class of size 1.
  label  = check_label(label, n)
  ulabel = unique(label)
  for (i in 1:length(ulabel)){
    if (sum(label==ulabel[i])==1){
      stop("* do.eslpp : no degerate class of size 1 is allowed.")
    }
  }
  if (any(is.na(label))||(any(is.infinite(label)))){warning("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.eslpp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. numk
  numk = as.integer(numk)
  if (!check_NumMM(numk,1,n/2,compact=FALSE)){stop("* do.eslpp : 'numk' should be an integer in [2,nrow(X)/2).")}
  #   5. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. preprocessing
  if (algpreprocess=="null"){
    trfinfo = list()
    trfinfo$type = "null"
    pX = as.matrix(X,nrow=nrow(X));
  } else {
    tmplist = aux.preprocess(X,type=algpreprocess)
    trfinfo = tmplist$info
    pX      = tmplist$pX
  }
  trfinfo$algtype = "linear"

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
  Stmp = array(0,c(n,n))
  for (i in 1:numk){
    tgtidx = clustidx[[i]]
    Stmp[tgtidx,tgtidx] = method_trfextlpp(PD[tgtidx,tgtidx],veca[i],vecb[i])
  }
  diag(Stmp) = 0.0

  ############# EXTENDED "SUPERVISED" SENSE
  S = array(0,c(n,n))
  for (i in 1:length(ulabel)){
    tgtidx = which(label==ulabel[i])
    S[tgtidx,tgtidx] = Stmp[tgtidx,tgtidx]
  }

  #   5. graph laplaciana and generalized eigenvalue problem
  D = diag(rowSums(S))
  L = D-S

  LHS = t(pX)%*%L%*%pX
  RHS = t(pX)%*%D%*%pX

  #   6. compute Projection Matrix
  geigs = geigen::geigen(LHS, RHS, TRUE)
  projection = matrix(geigs$vectors[,1:ndim],nrow=p)
  eigenvalue = as.vector(geigs$values[1:ndim])

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$eigval = eigenvalue
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}

