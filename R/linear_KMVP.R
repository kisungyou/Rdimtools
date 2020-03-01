#' Kernel-Weighted Maximum Variance Projection
#'
#' Kernel-Weighted Maximum Variance Projection (KMVP) is a generalization of
#' Maximum Variance Projection (MVP). Even though its name contains \emph{kernel}, it is
#' not related to kernel trick well known in the machine learning community. Rather, it
#' generalizes the binary penalization on class discrepancy,
#' \deqn{S_{ij} = \exp(-\|x_i-x_j\|^2/t) \quad\textrm{if}\quad C_i \ne C_j}
#' where \eqn{x_i} is an \eqn{i}-th data point and \eqn{t} a kernel bandwidth (\code{bandwidth}). \bold{Note} that
#' when the bandwidth value is too small, it might suffer from numerical instability and rank deficiency due to its formulation.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param bandwidth bandwidth parameter for heat kernel as the equation above.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \dontrun{
#' ## load iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.factor(iris$Species)
#'
#' ## perform KMVP with different bandwidths
#' out1 = do.kmvp(X, label, bandwidth=0.1)
#' out2 = do.kmvp(X, label, bandwidth=1)
#' out3 = do.kmvp(X, label, bandwidth=10)
#'
#' ## visualize
#' opar <- par(mfrow=c(1,3), no.readonly=TRUE)
#' plot(out1$Y, main="bandwidth=0.1", col=label)
#' plot(out2$Y, main="bandwidth=1",   col=label)
#' plot(out3$Y, main="bandwidth=10",  col=label)
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zhang_maximum_2007}{Rdimtools}
#'
#' @seealso \code{\link{do.mvp}}
#' @author Kisung You
#' @rdname linear_KMVP
#' @export
do.kmvp <- function(X, label, ndim=2,
                    preprocess=c("center","scale","cscale","decorrelate","whiten"),
                    bandwidth=1.0){
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
      stop("* do.kmvp : no degerate class of size 1 is allowed.")
    }
  }
  N = length(ulabel)
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.kmvp : 'ndim' is a positive integer in [1,#(covariates)).")}
  if (ndim>=N){
    stop("* do.kmvp : the method requires {ndim <= N-1}, where N is the number of classes.")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. bandwidth
  bandwidth = as.double(bandwidth)
  if (!check_NumMM(bandwidth,0,1e+10,compact=TRUE)){stop("* do.kmvp : 'bandwidth' should be a nonnegative real number.")}

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. preprocess of data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. perform PCA onto (N-1) dimensional space
  outPCA = do.pca(pX,ndim=(N-1))
  projection_first = aux.adjprojection(outPCA$projection)
  ppX = outPCA$Y

  #   3. Start of MVP algorithm here.
  #   3-1. Laplacian Part
  S = kmvp_S(ppX, label, bandwidth) # S_ij = 1 if Ci != Cj
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
    W[i,tgtidx] = kmvp_Gvec(dataw, tgtdata)
  }
  #   3-3. compute M
  Mhalf = diag(n)-W
  M     = (t(Mhalf)%*%Mhalf)

  #   4. solve geigen
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
kmvp_S <- function(X, label, bd){
  n = nrow(X)
  S = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (label[i]!=label[j]){
        diffvec = as.vector(X[i,]-X[j,])
        thevalue = exp(-sum(diffvec*diffvec)/bd)
        S[i,j] = thevalue
        S[j,i] = thevalue
      }
    }
  }
  return(S)
}
#' @keywords internal
#' @noRd
kmvp_Gvec <- function(vec, mat){
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
