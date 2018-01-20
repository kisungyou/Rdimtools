#' Maximum Margin Criterion
#'
#' Maximum Margin Criterion (MMC) is a linear supervised dimension reduction method that
#' maximizes average margin between classes. The cost function is defined as
#' \deqn{trace(S_b - S_w)}
#' where \eqn{S_b} is an overall variance of class mean vectors, and \eqn{S_w} refers to
#' spread of every class. Note that Principal Component Analysis (PCA) maximizes
#' total scatter, \eqn{S_t = S_b + S_w}.
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center" and two other options "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## generate 3 different groups of data X and label vector
#' x1 = matrix(rnorm(4*10), nrow=10)-20
#' x2 = matrix(rnorm(4*10), nrow=10)
#' x3 = matrix(rnorm(4*10), nrow=10)+20
#' X  = rbind(x1, x2, x3)
#' label = c(rep(1,10), rep(2,10), rep(3,10))
#'
#' ## perform MVP with different preprocessings
#' out1 = do.mmc(X, label, ndim=2)
#' out2 = do.mmc(X, label, ndim=2, preprocess="decorrelate")
#'
#' ## visualize
#' par(mfrow=c(1,2))
#' plot(out1$Y[,1], out1$Y[,2], main="centering")
#' plot(out2$Y[,1], out2$Y[,2], main="decorrelating")
#'
#' @references
#' \insertRef{li_efficient_2006}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_MMC
#' @export
do.mmc <- function(X, label, ndim=2, preprocess=c("center","decorrelate","whiten")){
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
      stop("* do.mmc : no degerate class of size 1 is allowed.")
    }
  }
  N = length(ulabel)
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.mmc : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocess of data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
  #   2. vector of proportion : pi
  nlabel = length(ulabel) # number of classes
  proportion = rep(0,nlabel)
  for (i in 1:nlabel){
    proportion[i] = (length(which(label==ulabel[i]))/n)
  }
  #   3. per-class and overall : mean vectors
  mean_PerClass = array(0,c(nlabel,p))
  for (i in 1:nlabel){
    idxlabel = which(label==ulabel[i])
    mean_PerClass[i,] = as.vector(colMeans(pX[idxlabel,]) )
  }
  mean_Overall = as.vector(colMeans(pX))
  #   4. per-class and overall : scatter
  scatter_PerClass = list()
  for (i in 1:nlabel){
    idxlabel = which(label==ulabel[i])
    scatter_PerClass[[i]] = mmc_scatter(pX[idxlabel,])
  }
  #   5. compute Sw
  Sw = array(0,c(p,p))
  for (i in 1:nlabel){
    Sw = Sw + ((proportion[i])*(scatter_PerClass[[i]]))
  }
  #   6. compute Sb
  Sb = array(0,c(p,p))
  for (i in 1:nlabel){
    vecdiff = as.vector(mean_PerClass[i,])-as.vector(mean_Overall)
    Sb = Sb + ((proportion[i])*outer(vecdiff,vecdiff))
  }
  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION
  costS = Sb-Sw
  projection = aux.adjprojection(RSpectra::eigs(costS, ndim)$vectors)

  #------------------------------------------------------------------------
  ## RETURN THE RESULTS
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}


#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
mmc_scatter <- function(tgtmat){
  N      = nrow(tgtmat)
  onesN  = rep(1,N)
  Cn     = diag(N)-(outer(onesN,onesN)/N)
  output = t(tgtmat)%*%Cn%*%tgtmat
  return(output)
}
