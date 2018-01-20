#' Kernel Maximum Margin Criterion
#'
#' Kernel Maximum Margin Criterion (KMMC) is a nonlinear variant of MMC method using kernel trick.
#' For computational simplicity, only the gaussian kernel is used with bandwidth parameter \code{t}.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center" and two other options "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param t bandwidth parameter for heat kernel in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
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
#' out1 = do.kmmc(X, label, t=0.1)
#' out2 = do.kmmc(X, label, t=1.0)
#' out3 = do.kmmc(X, label, t=10.0)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="bandwidth=0.1")
#' plot(out2$Y[,1], out2$Y[,2], main="bandwidth=1")
#' plot(out3$Y[,1], out3$Y[,2], main="bandwidth=10.0")
#'
#' @references
#' \insertRef{li_efficient_2006}{Rdimtools}
#'
#' @author Kisung You
#' @seealso \code{\link{do.mmc}}
#' @rdname nonlinear_KMMC
#' @export
do.kmmc <- function(X, label, ndim=2, preprocess=c("center","decorrelate","whiten"), t=1.0){
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
      stop("* do.kmmc : no degerate class of size 1 is allowed.")
    }
  }
  N = length(ulabel)
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.kmmc : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. t : kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t, 0, 1e+10, compact=FALSE)){stop("* do.kmmc : 't' is a bandwidth parameter for gaussian kernel.")}
  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocess of data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "nonlinear"
  #   2. vector of counts and proportion
  nlabel = length(ulabel) # number of classes
  vec_counts     = rep(0,nlabel)
  vec_proportion = rep(0,nlabel)
  for (i in 1:nlabel){
    vec_counts[i]     = length(which(label==ulabel[i]))
    vec_proportion[i] = (vec_counts[i])/n
  }
  #   3. Compute St
  onesn      = rep(1,n)
  kerneltype = c("gaussian",t)
  tmpKSt     = exp(-(as.matrix(dist(pX))^2)/(2*(t^2)))
  tilde_St   = (tmpKSt%*%(diag(n) - (outer(onesn,onesn)/n))%*%t(tmpKSt))/n
  #   4. m_i
  tilde_mi = array(0,c(nlabel,n))
  for (i in 1:nlabel){
    #   4-1. select class
    tgtidx  = which(label==ulabel[i])
    partmat = pX[tgtidx,]
    if (length(tgtidx)==1){
      partmat = matrix(partmat,nrow=1)
    }
    #   4-2. compute m_i with RCPP, YEAH!
    tilde_mi[i,] = as.vector(method_kmmcvec(pX, partmat, t))
  }
  #   4-3. Compute Sb
  tilde_Sb = array(0,c(n,n))
  mean_all = rep(0,n)
  for (i in 1:nlabel){
    mean_all = mean_all + ((vec_proportion[i])*as.vector(tilde_mi[i,]))
  }
  for (i in 1:nlabel){
    vecdiff = as.vector(tilde_mi[i,])-as.vector(mean_all)
    tilde_Sb = tilde_Sb + outer(vecdiff,vecdiff)*(vec_proportion[i])
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART
  costS = (2*tilde_Sb)-tilde_St
  eigsS = RSpectra::eigs(costS, ndim)

  #------------------------------------------------------------------------
  ## RETURN THE RESULTS
  result = list()
  result$Y = t(tmpKSt)%*%eigsS$vectors
  result$trfinfo = trfinfo
  return(result)
}



# -------------------------------------------------------------------------
# Even though this method is 'nonlinear', at first I misunderstood it as
# linear method in that its CPP function is listed in 'linear' graoup of Rcpp functions.

