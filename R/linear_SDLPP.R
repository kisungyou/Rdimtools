#' Sample-Dependent Locality Preserving Projection
#'
#' Many variants of Locality Preserving Projection are contingent on
#' graph construction schemes in that they sometimes return a range of
#' heterogeneous results when parameters are controlled to cover a wide range of values.
#' This algorithm takes an approach called \emph{sample-dependent construction} of
#' graph connectivity in that it tries to discover intrinsic structures of data
#' solely based on data.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param t kernel bandwidth in \eqn{(0,\infty)}.
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
#' @seealso \code{\link{do.lpp}}
#'
#' @examples
#' \donttest{
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.integer(iris$Species)
#'
#' ## compare with PCA
#' out1 <- do.pca(X,ndim=2)
#' out2 <- do.sdlpp(X, t=0.1)
#' out3 <- do.sdlpp(X, t=1)
#' out4 <- do.sdlpp(X, t=10)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2))
#' plot(out1$Y, col=label, main="PCA")
#' plot(out2$Y, col=label, main="SDLPP::t=0.1")
#' plot(out3$Y, col=label, main="SDLPP::t=1")
#' plot(out4$Y, col=label, main="SDLPP::t=10")
#' par(opar)
#' }
#'
#'
#' @references
#' \insertRef{yang_sampledependent_2010}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_SDLPP
#' @concept linear_methods 
#' @export
do.sdlpp <- function(X, ndim=2, t = 1.0,
                     preprocess=c("center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.sdlpp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. t
  t = as.double(t)
  if (!check_NumMM(t,0,1e+10,compact=FALSE)){stop("* do.sdlpp : 't' should be a positive real number.")}
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. preprocessing of data matrix
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. preliminary computations
  #   2-1. compute D : normalized squared distance
  D = (as.matrix(dist(pX, method="euclidean"))^2)
  for (i in 1:n){
    sumvecD = sum(D[i,])
    D[i,] = D[i,]/sumvecD
  }
  #   2-2. compute Ss : exponentiated D squared                  :: Problem with 'zerodiag'.
  Ss = exp(-D/(2*(t^2)))
  # if (diagzero){                                               :: follow the method directly.
  #   diag(Ss) = 0
  # }
  #   2-3. compute Ws : conditionally
  Ws = array(0,c(n,n))
  rowMeansSs = rowMeans(Ss)
  for (i in 1:n){
    for (j in 1:n){
      if (Ss[i,j] > (rowMeansSs[i])){
        Ws[i,j] = Ss[i,j]
      }
    }
  }
  #   2-4. compute Ds and Ls
  Wtilde = Ws+t(Ws)
  Ds = diag(rowSums(Ws))+diag(colSums(Ws))
  Ls = Ds-Wtilde

  #   3. main LPP part
  LHS = t(pX)%*%Ls%*%pX
  RHS = t(pX)%*%Ds%*%pX

  #   4. compute Projection Matrix
  projection = aux.geigen(LHS, RHS, ndim, maximal=FALSE)

  #------------------------------------------------------------------------
  ## RETURN
  #   1. adjust projection
  projection = aux.adjprojection(projection)
  #   2. return
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
