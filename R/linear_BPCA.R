#' Bayesian Principal Component Analysis
#'
#' Bayesian PCA (BPCA) is a further variant of PCA in that it imposes prior and encodes
#' basis selection mechanism. Even though the model is fully Bayesian, \code{do.bpca}
#' faithfully follows the original paper by Bishop in that it only returns the mode value
#' of posterior as an estimate, in conjunction with ARD-motivated prior as well as
#' consideration of variance to be estimated. Unlike PPCA, it uses full basis and returns
#' relative weight for each base in that the smaller \eqn{\alpha} value is, the more likely
#' corresponding column vector of \code{mp.W} to be selected as potential basis.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an option for preprocessing the data. Default is "center". See also \code{\link{aux.preprocess}}
#' for more details.
#' @param reltol stopping criterion for iterative update for EM algorithm.
#' @param maxiter maximum number of iterations allowed for EM algorithm.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)}  whose columns are principal components.}
#' \item{mp.itercount}{the number of iterations taken for EM algorithm to converge.}
#' \item{mp.sigma2}{estimated \eqn{\sigma^2} value via EM algorithm.}
#' \item{mp.alpha}{length-\code{ndim-1} vector of relative weight for each base in \code{mp.W}.}
#' \item{mp.W}{an \eqn{(ndim\times ndim-1)} matrix from EM update.}
#' }
#'
#'
#' @seealso \code{\link{do.pca}}, \code{\link{do.ppca}}
#' @author Kisung You
#' @references
#' \insertRef{bishop_bayesian_1999}{Rdimtools}
#'
#'
#' @examples
#' \dontrun{
#' ## generate swiss roll data
#' X = aux.gensamples()
#'
#' ## Compare PCA and BPCA
#' out1  <- do.pca(X, ndim=2, preprocess="center")
#' out2  <- do.bpca(X, ndim=2, preprocess="center")
#'
#' ## Visualize
#' par(mfrow=c(1,2))
#' plot(out1$Y[,1], out1$Y[,2], main="PCA")
#' plot(out2$Y[,1], out2$Y[,2], main="BPCA")
#' }
#'
#' @rdname linear_BPCA
#' @export
do.bpca <- function(X, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"),
                    reltol=1e-4, maxiter=123){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  # 1. data X
  aux.typecheck(X)
  # 2. ndim
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>=ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.bpca : 'ndim' is a positive integer in [1,#(covariates)).")
  }
  # 3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  # 4. reltol and maxiter : I trust users.. Never
  mepsil = .Machine$double.eps
  if (is.na(reltol)||is.infinite(reltol)||(reltol<mepsil)||(reltol>=1)||(!is.numeric(reltol))){
    stop("* do.bpca : 'reltol' should be in [machine epsilon,1)")
  }
  reltol = as.double(reltol)
  if ((!is.numeric(maxiter))||(maxiter<5)||(is.na(maxiter))||(is.infinite(maxiter))){
    stop("* do.bpca : 'maxiter' is a positive integer greater than or equal to 5.")
  }
  maxiter = as.integer(maxiter)


  #------------------------------------------------------------------------
  ## COMPUTATION
  #   1. Preprocessing the data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. Run using Rcpp
  rcppbpca = method_bpca(t(pX), reltol, maxiter);

  #   3. we select alpha with smallest values only.
  smallidx = order(rcppbpca$alpha)[1:ndim]
  #      in that we find projection as required
  mlsig2 = rcppbpca$sig2
  mlW    = rcppbpca$W[,smallidx]
  M = (t(mlW)%*%mlW)+(diag(ncol(mlW))*mlsig2)
  SOL = aux.bicgstab(M, t(mlW), verbose=FALSE)
  projection = t(SOL$x)


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y          = pX%*%projection
  result$trfinfo    = trfinfo
  result$projection = projection
  result$mp.itercount = rcppbpca$itercount # number of iterations
  result$mp.sigma2  = rcppbpca$sig2
  result$mp.alpha   = rcppbpca$alpha
  result$mp.W       = rcppbpca$W
  return(result)
}

