#' t-distributed Stochastic Neighbor Embedding
#'
#' \eqn{t}-distributed Stochastic Neighbor Embedding (t-SNE) is a variant of Stochastic Neighbor Embedding (SNE)
#' that mimicks patterns of probability distributinos over pairs of high-dimensional objects on low-dimesional
#' target embedding space by minimizing Kullback-Leibler divergence. While conventional SNE uses gaussian
#' distributions to measure similarity, t-SNE, as its name suggests, exploits a heavy-tailed Student t-distribution.
#' For \code{do.tsne}, we implemented a naive version of t-SNE.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param perplexity desired level of perplexity; ranging [5,50].
#' @param eta learning parameter.
#' @param maxiter maximum number of iterations.
#' @param jitter level of white noise added at the beginning.
#' @param jitterdecay decay parameter in (0,1). The closer to 0, the faster artificial noise decays.
#' @param momentum level of acceleration in learning.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param pca whether to use PCA as preliminary step; \code{TRUE} for using it, \code{FALSE} otherwise.
#' @param pcaratio proportion of variances explained in finding PCA preconditioning. See also \code{\link{do.pca}} for more details.
#' @param pcascale a logical; \code{FALSE} for using Covariance, \code{TRUE} for using Correlation matrix. See also \code{\link{do.pca}} for more details.
#' @param symmetric a logical; \code{FALSE} to solve it naively, and \code{TRUE} to adopt symmetrization scheme.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate swiss roll data
#' X = aux.gensamples(n=100)
#'
#' ## compare different perplexity
#' out1 <- do.tsne(X, ndim=2, perplexity=5)
#' out2 <- do.tsne(X, ndim=2, perplexity=10)
#' out3 <- do.tsne(X, ndim=2, perplexity=50)
#'
#' ## Visualize three different projections
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="tSNE::perplexity=5")
#' plot(out1$Y[,1], out1$Y[,2], main="tSNE::perplexity=10")
#' plot(out1$Y[,1], out1$Y[,2], main="tSNE::perplexity=50")
#' }
#'
#' @seealso \code{\link{do.sne}}
#' @references
#' \insertRef{van_der_maaten_visualizing_2008}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_TSNE
#' @export
do.tsne <- function(X,ndim=2,perplexity=30,eta=0.05,maxiter=2000,
                    jitter=0.3,jitterdecay=0.99,momentum=0.5,
                    preprocess=c("null","center","scale","cscale","decorrelate","whiten"),
                    pca=TRUE,pcaratio=0.90,pcascale=FALSE,symmetric=FALSE){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  #   1-1. (integer) ndim
  if (!is.numeric(ndim)||(ndim<1)||(ndim>ncol(X))){
    stop("* do.tsne : 'ndim' is an integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  #   1-2. perplexity
  if (!is.numeric(perplexity)||is.na(perplexity)||is.infinite(perplexity)||(perplexity<=0)){
    stop("* do.tsne : perplexity should be a positive real number.")
  }
  if ((perplexity < 5)||(perplexity > 50)){
    message("* do.tsne : a desired perplexity value is in [5,50].")
  }

  # obsolete params.
  BarnesHut=FALSE
  BHtheta=0.5
  # 2. Input Parameters
  #   2-1. (double) eta = 0.5; learning parameter
  if (!is.numeric(eta)||is.na(eta)||is.infinite(eta)||(eta<=0)){
    stop("* do.tsne : learning rate 'eta' should be a positive real number.")
  }

  #   2-2. (integer) maxiter = 2000; maximum number of iterations
  if (!is.numeric(maxiter)||(maxiter<2)||(is.na(maxiter))||(is.infinite(maxiter))){
    stop("* do.tsne : maxiter should be suited for the number of iterations.")
  }
  #   2-3. (double) jitter = 0.3; random errors
  if (!is.numeric(jitter)||(is.na(jitter))||(is.infinite(jitter))||(jitter<0)){
    stop("* do.tsne : 'jitter' should be a positive real number.")
  }
  #   2-4. (double) jitterdecay = 0.99; decaying factor of jitter
  decay = jitterdecay
  if (!is.numeric(decay)||(is.na(decay))||(is.infinite(decay))||(decay<=0)||(decay>=1)){
    stop("* do.tsne : 'jitterdecay' is a real number between (0,1).")
  }
  #   2-5. (double) momentum = 0.5
  if ((!is.numeric(momentum))||(is.na(momentum))||(is.infinite(momentum))||(momentum<=0)){
    stop("* do.tsne : 'momentum' should be a positive real number.")
  }
  #   2-6. (char) preprocess = 'center'
  algpreprocess = match.arg(preprocess)
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2-7. (bool) pca = TRUE/FALSE
  #     If pca = TRUE
  #        pcaratio (0,1) : variance ratio
  #        pcascale       : TRUE/FALSE
  pcaflag = pca; if(!is.logical(pcaflag)){stop("* do.tsne : 'pca' is a logical variable.")}
  if (!is.numeric(pcaratio)||(pcaratio<=0)||(pcaratio>=1)||is.na(pcaratio)){
    stop("* do.tsne : pcaratio should be in (0,1).")
  }
  scaleflag = pcascale; if (!is.logical(scaleflag)){
    stop("* do.tsne : pcascale is either TRUE or FALSE.")
  }
  if (pcaflag){
    pcaout = do.pca(pX,ndim="auto",cor=scaleflag,preprocess="center",varratio=pcaratio)
    if (ncol(pcaout$Y)<=ndim){
      message("* do.tsne : PCA scaling has gone too far.")
      message("* do.tsne : Pass non-scaled data to t-SNE algortihm.")
      tpX = t(pX)
    } else {
      tpX = t(pcaout$Y)
    }
  } else {
    tpX = t(pX)
  }

  #   2-8. (bool) BarnesHut : TRUE/FALSE
  #               BHtheta   : 0 means exact tSNE that runs with my own code.
  BHflag = BarnesHut;
  if (!is.logical(BHflag)){
    stop("* do.tsne : 'BarnesHut' is a logical variable.")
  }
  if (!is.numeric(BHtheta)||is.na(BHtheta)||(BHtheta<0)||is.infinite(BHtheta)){
    stop("* do.tsne : BHtheta is invalid. It should be >= 0.")
  }
  BHtheta = as.double(BHtheta)

  # 3. Run Main Algorithm
  # if (!BHflag){
  #   Perp = aux_perplexity(tpX,perplexity);
  #   P = as.matrix(Perp$P)
  #   vars = as.vector(Perp$vars)
  #
  #   Y = t(as.matrix(method_tsne(P,ndim,eta,maxiter,jitter,decay,momentum)))
  # } else {
  #   pX = t(tpX)
  #   out = Rtsne(pX,dims=ndim,theta=BHtheta,perplexity=perplexity,pca=TRUE,max_iter=maxiter,
  #               momentum=momentum,eta=eta)
  #   Y = (out$Y)
  # }
  Perp = aux_perplexity(tpX,perplexity);
  P = as.matrix(Perp$P)
  vars = as.vector(Perp$vars)
  Y = t(as.matrix(method_tsne(P,ndim,eta,maxiter,jitter,decay,momentum)))

  # 5. result
  if (any(is.infinite(Y))||any(is.na(Y))){
    message("* do.tsne : t-SNE not successful; having either Inf or NA values.")
    result = NA
    return(result)
  } else {
    result = list()
    result$Y = Y
    result$trfinfo = trfinfo
  }
  return(result)
}
