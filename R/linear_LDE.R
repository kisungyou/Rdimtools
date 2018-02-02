#' Local Discriminant Embedding
#'
#' Local Discriminant Embedding (LDE) is a supervised algorithm that learns
#' the embedding for the submanifold of each class. Its idea is to same-class
#' data points maintain their original neighborhood information while
#' segregating different-class data distinct from each other.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param t kernel bandwidth in \eqn{(0,\infty)}.
#' @param numk the number of neighboring points for k-nn graph construction.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center" and other options of "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @references
#' \insertRef{hwann-tzong_chen_local_2005}{Rdimtools}
#'
#' @examples
#' ## generate data of 2 types with clear difference
#' diff = 15
#' dt1  = aux.gensamples(n=123)-diff;
#' dt2  = aux.gensamples(n=123)+diff;
#'
#' ## merge the data and create a label correspondingly
#' Y      = rbind(dt1,dt2)
#' label  = c(rep(1,123), rep(2,123))
#'
#' ## try different neighborhood size
#' out1 <- do.lde(Y, label, numk=5)
#' out2 <- do.lde(Y, label, numk=10)
#' out3 <- do.lde(Y, label, numk=25)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="LDE::k=5")
#' plot(out2$Y[,1], out2$Y[,2], main="LDE::k=10")
#' plot(out3$Y[,1], out3$Y[,2], main="LDE::k=25")
#'
#' @author Kisung You
#' @rdname linear_LDE
#' @export
do.lde <- function(X, label, ndim=2, t=1.0, numk=max(ceiling(nrow(X)/10),2),
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
      stop("* do.lde : no degerate class of size 1 is allowed.")
    }
  }
  if (any(is.na(label))||(any(is.infinite(label)))){warning("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.lde : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. t
  t = as.double(t)
  if (!check_NumMM(t,.Machine$double.eps,Inf,compact=TRUE)){stop("* do.lde : 't' should be a positive real number.")}
  #   5. numk
  numk = as.integer(numk)
  if (!check_NumMM(numk,1,n/2,compact=FALSE)){stop("* do.lde : 'numk' should be an integer in [2,nrow(X)/2).")}
  #   6. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   Pre. preprocessing of data matrix
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"

  #   1. construct G1 (original G) and G2 (G') for same-class and different-class connectivty
  #   1-1. find k-neighborhood graph
  nbdtype   = c("knn",numk)
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric="union")
  Dmask     = nbdstruct$mask
  #   1-2. logical based-on class information
  conn_same = lde_perclass_logical(label)
  conn_diff = 1-conn_same
  #   1-3. connectivity
  G1 = Dmask*conn_same
  G2 = Dmask*conn_diff

  #   2. build AFFINITY matrix
  expD = exp(-(as.matrix(dist(pX))^2)/t)
  W1 = expD*G1
  W2 = expD*G2

  #   3. Want To Find Embedding
  LHS = t(pX)%*%(diag(rowSums(W2))-W2)%*%pX
  RHS = t(pX)%*%(diag(rowSums(W1))-W1)%*%pX

  #   4. compute Projection Matrix : use lowest ones
  projection = aux.geigen(LHS, RHS, ndim, maximal=FALSE)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}



#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
lde_perclass_logical <- function(label){
  n = length(label)
  out1 = matrix(0,nrow=n,ncol=n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (label[i]==label[j]){
        out1[i,j] = 1
        out1[j,i] = 1
      }
    }
  }
  diag(out1) = 1
  return(out1)
}
