#' Kernel-Weighted Unsupervised Discriminant Projection
#'
#' Kernel-Weighted Unsupervised Discriminant Projection (KUDP) is a generalization of UDP where
#' proximity is given by weighted values via heat kernel,
#' \deqn{K_{i,j} = \exp(-\|x_i-x_j\|^2/kbandwidth)}
#' whence UDP uses binary connectivity. If \code{kbandwidth} is \eqn{+\infty}, it becomes
#' a standard UDP problem. Like UDP, it also performs PCA preprocessing for rank-deficient case.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center" and two other options "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param kbandwidth bandwidth parameter for heat kernel as the equation above.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' \item{interimdim}{the number of PCA target dimension used in preprocessing.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate sample data
#' X = aux.gensamples(n=200)
#'
#' ## use different kernel bandwidth
#' out1 <- do.kudp(X, kbandwidth=0.1)
#' out2 <- do.kudp(X, kbandwidth=10)
#' out3 <- do.kudp(X, kbandwidth=1000)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1],out1$Y[,2],main="bandwidth=0.1")
#' plot(out2$Y[,1],out2$Y[,2],main="bandwidth=10")
#' plot(out3$Y[,1],out3$Y[,2],main="bandwidth=1000")
#' }
#'
#' @seealso \code{\link{do.udp}}
#' @author Kisung You
#' @rdname linear_KUDP
#' @references
#' \insertRef{yang_globally_2007}{Rdimtools}
#'
#' @export
do.kudp <- function(X, ndim=2, type=c("proportion",0.1),
                   preprocess=c("center","decorrelate","whiten"),
                   kbandwidth=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  M = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.kudp : 'ndim' is a positive integer in [1,#(covariates)].")}
  #   3. preprocessing
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  trfinfo$algtype = "linear"
  pX      = tmplist$pX
  #   4. kernel-weights
  kbandwidth = as.double(kbandwidth)
  if (!check_NumMM(kbandwidth,0,1e+10,compact=TRUE)){stop("* do.kudp : 'kbandwidth' should be a nonnegative real number.")}

  #   4. neighborhood creation
  nbdtype = type
  nbdsymmetric = "intersect"
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  #   5. kernelized
  K = exp(-(as.matrix(dist(pX))^2)/kbandwidth); diag(K) = 0; ## this part i'm not sure
  H = K*(nbdstruct$mask)
  D = diag(colSums(H))
  L = D-H                  # graph laplacian / local scatter kernel

  #------------------------------------------------------------------------
  ## INTERMEDIATE PCA
  # 1. compute St
  tmpSt = kudp_ST(pX)
  # 2. target rank
  tmpndim = min((max(as.integer(Matrix::rankMatrix(tmpSt)), (ndim+1))), p)
  # 3. perform PCA
  if (tmpndim==p){
    proj_first = diag(rep(1,p))
    tmp_X = pX

    ## MAIN PART for No Need To Suffer from Low-Dimensional Issue
    # 1. compute Sn, Sl from St
    final_ST = kudp_ST(tmp_X)
    final_SL = (t(tmp_X)%*%L%*%(tmp_X))/(M*M)
    final_SN = final_ST - final_SL

    # 2. geigen : use bottom ones
    proj_second = aux.geigen(final_SN, final_SL, ndim, maximal=TRUE)
    proj_all = (proj_first %*% proj_second)
  } else {
    eigSt = eigen(tmpSt)
    topeigvals = eigSt$values[1:tmpndim]
    proj_first = eigSt$vectors[,1:tmpndim] # P : (p-by-tmpndim)
    Xtilde = pX%*%proj_first

    ## MAIN PART for Low-Dimensional Case : $3.3. UDP Algorithm
    # Step 3.
    ST_tilde = diag(topeigvals)
    # Step 4.
    SL_tilde = (t(Xtilde)%*%L%*%Xtilde)/(M*M)
    SN_tilde = (ST_tilde-SL_tilde)

    proj_second = aux.geigen(SN_tilde, SL_tilde, ndim, maximal=TRUE)
    proj_all    = (proj_first%*%proj_second)
  }

  #------------------------------------------------------------------------
  ## RETURN THE RESULTS

  # 2. result list
  result = list()
  result$Y = pX%*%proj_all
  result$trfinfo = trfinfo
  result$projection = proj_all
  result$interimdim = tmpndim

  # 3. return
  return(result)
}



# auxiliary functions for UDP ---------------------------------------------
#' @keywords internal
#' @noRd
kudp_ST <- function(X){
  M = nrow(X)
  p = ncol(X)
  output = array(0,c(p,p))
  for (i in 1:(M-1)){
    veci = X[i,]
    for (j in (i+1):M){
      vecj = X[j,]
      vecdiff = (veci-vecj)
      output = output + outer(vecdiff, vecdiff)
    }
  }
  output = output/(M*M)
  return(output)
}
