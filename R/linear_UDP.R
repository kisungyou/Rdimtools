#' Unsupervised Discriminant Projection
#'
#' Unsupervised Discriminant Projection (UDP) aims finding projection that balances local and global scatter.
#' Even though the name contains the word \emph{Discriminant}, this algorithm is \emph{unsupervised}. The
#' term there reflects its algorithmic tactic to discriminate distance points not in the neighborhood of each data point.
#' It performs PCA as intermittent preprocessing for rank singularity issue. Authors clearly mentioned that it is inspired
#' by Locality Preserving Projection, which minimizes the local scatter only.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
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
#' ## load iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' ## use different connectivity level
#' out1 <- do.udp(X, type=c("proportion",0.05))
#' out2 <- do.udp(X, type=c("proportion",0.10))
#' out3 <- do.udp(X, type=c("proportion",0.25))
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, pch=19, main="connectivity 5%")
#' plot(out2$Y, col=label, pch=19, main="connectivity 10%")
#' plot(out3$Y, col=label, pch=19, main="connectivity 25%")
#' par(opar)
#'
#' @author Kisung You
#' @rdname linear_UDP
#' @references
#' \insertRef{yang_globally_2007}{Rdimtools}
#'
#' @seealso \code{\link{do.lpp}}
#' @concept linear_methods
#' @export
do.udp <- function(X, ndim=2, type=c("proportion",0.1), preprocess=c("center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  M = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.udp : 'ndim' is a positive integer in [1,#(covariates)].")}
  #   3. preprocessing
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #   4. neighborhood creation
  nbdtype = type
  nbdsymmetric = "intersect"
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  H = (nbdstruct$mask)*1.0
  D = diag(colSums(H))
  L = D-H                  # graph laplacian / local scatter kernel

  #------------------------------------------------------------------------
  ## INTERMEDIATE PCA
  # 1. compute St
  tmpSt = udp_ST(pX)
  # 2. target rank
  tmpndim = min((max(round(aux_rank(tmpSt)), (ndim+1))), p)
  # 3. perform PCA
  if (tmpndim==p){
    proj_first = diag(rep(1,p))
    tmp_X = pX

    ## MAIN PART for No Need To Suffer from Low-Dimensional Issue
    # 1. compute Sn, Sl from St
    final_ST = udp_ST(tmp_X)
    final_SL = (t(tmp_X)%*%L%*%(tmp_X))/(M*M)
    final_SN = final_ST - final_SL

    # 2. gEVD : ascending
    #    now we don't need this part.

    # 3. new projection : find the largest/maximal ones
    proj_second = aux.geigen(final_SN, final_SL, ndim, maximal=TRUE)
    proj_all = (proj_first %*% proj_second)
  } else {
    eigSt = eigen(tmpSt)
    topeigvals = eigSt$values[1:tmpndim]
    proj_first = aux.adjprojection(eigSt$vectors[,1:tmpndim]) # P : (p-by-tmpndim)
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
  # 1. adjust projection
  proj_all = aux.adjprojection(proj_all)

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
udp_ST <- function(X){
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
