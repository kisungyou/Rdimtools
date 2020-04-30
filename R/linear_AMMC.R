#' Adaptive Maximum Margin Criterion
#'
#' Adaptive Maximum Margin Criterion (AMMC) is a supervised linear dimension reduction method.
#' The method uses different weights to characterize the different contributions of the
#' training samples embedded in MMC framework. With the choice of  \code{a=0}, \code{b=0}, and
#' \code{lambda=1}, it is identical to standard MMC method.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param a tuning parameter for between-class weight in \eqn{[0,\infty)}.
#' @param b tuning parameter for within-class weight in \eqn{[0,\infty)}.
#' @param lambda balance parameter for between-class and within-class scatter matrices in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
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
#' ## try different lambda values
#' out1 = do.ammc(X, label, lambda=0.1)
#' out2 = do.ammc(X, label, lambda=1)
#' out3 = do.ammc(X, label, lambda=10)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="AMMC::lambda=0.1", pch=19, cex=0.5, col=label)
#' plot(out2$Y, main="AMMC::lambda=1",   pch=19, cex=0.5, col=label)
#' plot(out3$Y, main="AMMC::lambda=10",  pch=19, cex=0.5, col=label)
#' par(opar)
#'
#' @references
#' \insertRef{lu_adaptive_2011}{Rdimtools}
#'
#' @seealso \code{\link{do.mmc}}
#' @author Kisung You
#' @rdname linear_AMMC
#' @concept linear_methods
#' @export
do.ammc <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"),
                    a=1.0, b=1.0, lambda=1.0){
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
      stop("* do.ammc : no degerate class of size 1 is allowed.")
    }
  }
  nlabel = length(ulabel)
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.ammc : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. a, b : tuning parameters
  a = as.double(a)
  b = as.double(b)
  lambda = as.double(lambda)
  if (!check_NumMM(a,0,1e+10,compact=TRUE)){stop("* do.ammc : 'a' should be a nonnegative real number.")}
  if (!check_NumMM(b,0,1e+10,compact=TRUE)){stop("* do.ammc : 'b' should be a nonnegative real number.")}
  if (!check_NumMM(lambda,0,Inf,compact=FALSE)){stop("* do.ammc : 'lambda' should be a positive real number.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocess of data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. per-class and overall : mean vectors
  meanvectors   = ammc_meanvec(pX, label, ulabel)
  mean_Overall  = meanvectors$overall
  mean_PerClass = meanvectors$class
  #   3. adaptive scatter matrices
  adaSb = ammc_adaSb(mean_PerClass, a)
  adaSw = ammc_adaSw(pX, label, b)

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION
  costS = (adaSb - (lambda*adaSw))
  projection = aux.adjprojection(RSpectra::eigs(costS, ndim)$vectors)

  #------------------------------------------------------------------------
  ## RETURN THE RESULTS
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}



# auxiliary for AMMC ------------------------------------------------------
#' @keywords internal
#' @noRd
ammc_meanvec <- function(X, label, ulabel){
  p      = ncol(X)
  nlabel = length(ulabel)

  mean_Overall = as.vector(colMeans(X))
  mean_Class   = array(0,c(nlabel,p))
  for (i in 1:nlabel){
    idxlabel = which(label==ulabel[i])
    mean_Class[i,] = as.vector(colMeans(X[idxlabel,]))
  }

  output = list()
  output$overall = mean_Overall
  output$class   = mean_Class
  return(output)
}
#' @keywords internal
#' @noRd
ammc_adaSb <- function(mat, a){
  c  = nrow(mat)
  p  = ncol(mat)
  Sb = array(0,c(p,p))
  for (i in 1:c){
    vec1 = as.vector(mat[i,])
    for (j in 1:c){
      vec2 = as.vector(mat[j,])
      if (j!=i){
        vecdiff = (vec1-vec2)
        weight  = ((sqrt(sum(vecdiff*vecdiff)))^(-a))
        Sb = Sb + weight*outer(vecdiff,vecdiff)
      }
    }
  }
  return(Sb)
}
#' @keywords internal
#' @noRd
ammc_adaSw <- function(X, label, b){
  n = nrow(X)
  p = ncol(X)
  if (length(label)!=n){
    stop("* ammc_adaSw.")
  }
  ulabel = unique(label)
  c      = length(ulabel)
  Sw = array(0,c(p,p))
  for (i in 1:c){
    idxlabel = which(label==ulabel[i])
    ni = length(idxlabel)
    mi = as.vector(colMeans(X[idxlabel,]))
    for (j in 1:ni){
      cidx = as.integer(idxlabel[j])
      cvec = as.vector(X[cidx,])

      vecdiff = cvec-mi
      weight  = ((sqrt(sum(vecdiff*vecdiff)))^b)
      Sw      = Sw + weight*outer(vecdiff,vecdiff)
    }
  }
  return(Sw)
}
