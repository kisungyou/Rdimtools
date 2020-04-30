#' Maximum Variance Projection
#'
#' Maximum Variance Projection (MVP) is a supervised method based on linear discriminant analysis (LDA).
#' In addition to classical LDA, it further aims at preserving local information by capturing
#' the local geometry of the manifold via the following proximity coding,
#' \deqn{S_{ij} = 1\quad\textrm{if}\quad C_i \ne C_j\quad\textrm{and} = 0 \quad\textrm{otherwise}},
#' where \eqn{C_i} is the label of an \eqn{i}-th data point.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
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
#' \donttest{
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150, 50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' ## perform MVP with different preprocessings
#' out1 = do.mvp(X, label)
#' out2 = do.mvp(X, label, preprocess="decorrelate")
#' out3 = do.mvp(X, label, preprocess="whiten")
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, pch=19, main="centering")
#' plot(out2$Y, col=label, pch=19, main="decorrelating")
#' plot(out3$Y, col=label, pch=19, main="whitening")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zhang_maximum_2007}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_MVP
#' @concept linear_methods
#' @export
do.mvp <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten")){
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
      stop("* do.mvp : no degerate class of size 1 is allowed.")
    }
  }
  N = length(ulabel)
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.mvp : 'ndim' is a positive integer in [1,#(covariates)).")}
  if (ndim>=N){
    stop("* do.mvp : the method requires {ndim <= N-1}, where N is the number of classes.")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. preprocess of data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. perform PCA onto (N-1) dimensional space
  outPCA = do.pca(pX,ndim=(N-1))
  projection_first = outPCA$projection
  ppX = outPCA$Y

  #   3. Start of MVP algorithm here.
  #   3-1. Laplacian Part
  S = mvp_S(ppX, label) # S_ij = 1 if Ci != Cj
  D = diag(rowSums(S))
  L = D-S
  #   3-2. LLE approximation part # we need to construct W !!
  W = array(0,c(n,n))
  for (i in 1:n){
    # 1. the data
    dataw = ppX[i,]
    # 2. target index and corresponding data
    #    note that the index 'i' should be removed
    tgtidx  = setdiff(which(label==label[i]),i)
    tgtdata = ppX[tgtidx,]
    # 3. compute w and assign
    W[i,tgtidx] = mvp_Gvec(dataw, tgtdata)
  }
  #   3-3. compute M
  Mhalf = diag(n)-W
  M     = (t(Mhalf)%*%Mhalf)

  #   4. solve geigen with lowest
  LHS = t(ppX)%*%M%*%ppX
  RHS = t(ppX)%*%L%*%ppX
  projection_second = aux.geigen(LHS,RHS,ndim,maximal=FALSE)
  projection_all = projection_first%*%projection_second

  #------------------------------------------------------------------------
  ## RETURN THE RESULTS
  result = list()
  result$Y = pX%*%projection_all
  result$trfinfo = trfinfo
  result$projection = projection_all
  return(result)
}



#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
mvp_S <- function(X, label){
  n = nrow(X)
  S = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (label[i]!=label[j]){
        S[i,j] = 1
        S[j,i] = 1
      }
    }
  }
  return(S)
}
#' @keywords internal
#' @noRd
mvp_Gvec <- function(vec, mat){
  # 1. settings
  n = nrow(mat)
  G = array(0,c(n,n))
  # 2. iterate to construct
  for (i in 1:n){
    veci = mat[i,]
    for (j in i:n){
      vecj = mat[j,]

      diffi = vec-veci
      diffj = vec-vecj
      if (i==j){
        G[i,i] = sum(diffi*diffj)
      } else {
        valueout = sum(diffi*diffj)
        G[i,j] = valueout
        G[j,i] = valueout
      }
    }
  }
  # 3. compute the vector
  Ginv = 1/G
  if (any(is.infinite(Ginv))){
    Ginv[which(is.infinite(Ginv))] = 0
  }
  denominator = sum(Ginv)
  woutput = (rowSums(Ginv))/(sum(Ginv))
  return(woutput)
}
