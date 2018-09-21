#' Stochastic Neighbor Embedding
#'
#' Stochastic Neighbor Embedding (SNE) is a probabilistic approach to mimick distributional
#' description in high-dimensional - possible, nonlinear - subspace on low-dimensional target space.
#' \code{do.sne} fully adopts algorithm details in an original paper by Hinton and Roweis (2002).
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param perplexity desired level of perplexity; ranging [5,50].
#' @param eta learning parameter.
#' @param maxiter maximum number of iterations.
#' @param jitter level of white noise added at the beginning.
#' @param jitterdecay decay parameter in \eqn{(0,1)}. The closer to 0, the faster artificial noise decays.
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
#' \item{vars}{a vector containing betas used in perplexity matching.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate ribbon-shaped data
#' ## in order to pass CRAN pretest, n is set to be small.
#' X = aux.gensamples(dname="ribbon",n=99)
#'
#' ## 1. pca scaling with 90% variances explained
#' output1 <- do.sne(X,ndim=2,pca=TRUE)
#'
#' ## 2. pca scaling wiht 50% variances explained
#' output2 <- do.sne(X,ndim=2,pca=TRUE,pcaratio=0.50)
#'
#' ## 3. Setting 2 + smaller level of perplexity
#' output3 <- do.sne(X,ndim=2,pca=TRUE,pcaratio=0.50,perplexity=10)
#'
#' ## Visualize three different projections
#' par(mfrow=c(1,3))
#' if ((length(output1)!=1)&&(!is.na(output1))){plot(output1$Y[,1],output1$Y[,2],main="Setting 1")}
#' if ((length(output1)!=1)&&(!is.na(output2))){plot(output2$Y[,1],output2$Y[,2],main="Setting 2")}
#' if ((length(output1)!=1)&&(!is.na(output3))){plot(output3$Y[,1],output3$Y[,2],main="Setting 3")}
#' }
#'
#' @references
#' \insertRef{hinton_stochastic_2003}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_SNE
#' @export
do.sne <- function(X,ndim=2,perplexity=30,eta=0.05,maxiter=2000,
                   jitter=0.3,jitterdecay=0.99,momentum=0.5,
                   preprocess=c("null","center","scale","scale","decorrelate","whiten"),
                   pca=TRUE,pcaratio=0.90,pcascale=FALSE,symmetric=FALSE){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)
  #   1-1. (integer) ndim
  if (!is.numeric(ndim)||(ndim<1)||(ndim>ncol(X))){
    stop("* do.sne : 'ndim' is an integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  #   1-2. perplexity
  if (!is.numeric(perplexity)||is.na(perplexity)||is.infinite(perplexity)||(perplexity<=0)){
    stop("* do.sne : perplexity should be a positive real number.")
  }
  if ((perplexity < 5)||(perplexity > 50)){
    message("* do.sne : a desired perplexity value is in [5,50].")
  }

  # 2. Input Parameters
  #   2-1. (double) eta = 0.5; learning parameter
  if (!is.numeric(eta)||is.na(eta)||is.infinite(eta)||(eta<=0)){
    stop("* do.sne : learning rate 'eta' should be a positive real number.")
  }

  #   2-2. (integer) maxiter = 2000; maximum number of iterations
  if (!is.numeric(maxiter)||(maxiter<2)||(is.na(maxiter))||(is.infinite(maxiter))){
    stop("* do.sne : maxiter should be suited for the number of iterations.")
  }
  #   2-3. (double) jitter = 0.3; random errors
  if (!is.numeric(jitter)||(is.na(jitter))||(is.infinite(jitter))||(jitter<0)){
    stop("* do.sne : 'jitter' should be a positive real number.")
  }
  #   2-4. (double) jitterdecay = 0.99; decaying factor of jitter
  decay = jitterdecay
  if (!is.numeric(decay)||(is.na(decay))||(is.infinite(decay))||(decay<=0)||(decay>=1)){
    stop("* do.sne : 'jitterdecay' is a real number between (0,1).")
  }
  #   2-5. (double) momentum = 0.5
  if ((!is.numeric(momentum))||(is.na(momentum))||(is.infinite(momentum))||(momentum<=0)){
    stop("* do.sne : 'momentum' should be a positive real number.")
  }

  algpreprocess = match.arg(preprocess)
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2-7. (bool) pca = TRUE/FALSE
  #     If pca = TRUE
  #        pcaratio (0,1) : variance ratio
  #        pcascale       : TRUE/FALSE
  pcaflag = pca; if(!is.logical(pcaflag)){stop("* do.sne : 'pca' is a logical variable.")}
  if (!is.numeric(pcaratio)||(pcaratio<=0)||(pcaratio>=1)||is.na(pcaratio)){
    stop("* do.sne : pcaratio should be in (0,1).")
  }
  scaleflag = pcascale; if (!is.logical(scaleflag)){
    stop("* do.sne : pcascale is either TRUE or FALSE.")
  }
  if (pcaflag){
    pcaout = do.pca(pX,ndim="auto",cor=scaleflag,preprocess="center",varratio=pcaratio)
    if (ncol(pcaout$Y)<=ndim){
      message("* do.sne : PCA scaling has gone too far.")
      message("* do.sne : Pass non-scaled data to SNE algortihm.")
      tpX = t(pX)
    } else {
      tpX = t(pcaout$Y)
    }
  } else {
    tpX = t(pX)
  }


  # 3. Compute Perplexity Matrix P from original data
  Perp = aux_perplexity(tpX,perplexity);
  P = as.matrix(Perp$P)
  vars = as.vector(Perp$vars)


  # 4. run main SNE algorithm : Symmetric Key is now used.
  if (symmetric==FALSE){
    Y = t(as.matrix(method_sne(P,ndim,eta,maxiter,jitter,decay,momentum)))
  } else if (symmetric==TRUE){
    Y = t(as.matrix(method_snesym(P,ndim,eta,maxiter,jitter,decay,momentum)))
  }

  # 5. result
  if (any(is.infinite(Y))||any(is.na(Y))){
    message("* do.tsne : t-SNE not successful; having either Inf or NA values.")
    result = NA
    return(result)
  } else {
    result = list()
    result$Y = Y
    result$trfinfo = trfinfo
    result$vars    = vars
  }
  return(result)
}
