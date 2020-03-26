#' Semi-Supervised Adaptive Maximum Margin Criterion
#'
#' Semi-Supervised Adaptive Maximum Margin Criterion (SAMMC) is a semi-supervised variant of
#' AMMC by making use of both labeled and unlabeled data.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param a tuning parameter for between-class weight in \eqn{[0,\infty)}.
#' @param b tuning parameter for within-class weight in \eqn{[0,\infty)}.
#' @param lambda balance parameter for between-class and within-class scatter matrices in \eqn{(0,\infty)}.
#' @param beta balance parameter for within-class scatter of the labeled data and consistency of the whole data in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \donttest{
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=33)-50
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+50
#'
#' ## merge the data and create a label correspondingly
#' X      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## copy a label and let 20% of elements be missing
#' nlabel = length(label)
#' nmissing = round(nlabel*0.20)
#' label_missing = label
#' label_missing[sample(1:nlabel, nmissing)]=NA
#'
#' ## try different balancing
#' out1 = do.sammc(X, label_missing, beta=0.1)
#' out2 = do.sammc(X, label_missing, beta=1)
#' out3 = do.sammc(X, label_missing, beta=10)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="SAMMC::beta=0.1")
#' plot(out2$Y, main="SAMMC::beta=1")
#' plot(out3$Y, main="SAMMC::beta=10")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{lu_adaptive_2011}{Rdimtools}
#'
#' @seealso \code{\link{do.mmc}}, \code{\link{do.ammc}}
#' @author Kisung You
#' @rdname linear_SAMMC
#' @export
do.sammc  <- function(X, label, ndim=2, type=c("proportion",0.1),
                      preprocess=c("center","scale","cscale","decorrelate","whiten"),
                      a=1.0, b=1.0, lambda=1.0, beta=1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label : check and return a de-factored vector
  label  = check_label(label, n)
  ulabel = unique(label)
  if (all(!is.na(ulabel))){
    message("* Semi-Supervised Learning : there is no missing labels. Consider using Supervised methods.")
  }
  if (any(is.infinite(ulabel))){
    stop("* Semi-Supervised Learning : Inf is not allowed in label.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.sammc : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. type
  nbdtype = type
  nbdsymmetric = "union"
  #   5. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   6. a and b : tuning parameters
  a = as.double(a)
  b = as.double(b)
  if (!check_NumMM(a,0,1e+10,compact=TRUE)){stop("* do.sammc : 'a' should be a nonnegative real number.")}
  if (!check_NumMM(b,0,1e+10,compact=TRUE)){stop("* do.sammc : 'b' should be a nonnegative real number.")}
  #   7. lambda and beta
  lambda = as.double(lambda)
  beta   = as.double(beta)
  if (!check_NumMM(lambda,0,Inf,compact=FALSE)){stop("* do.sammc : 'lambda' should be a positive real number.")}
  if (!check_NumMM(beta,0,Inf,compact=FALSE)){stop("* do.sammc : 'beta' should be a positive real number.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocess of data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #   2. non-missing labels
  idxfilled   = which(!is.na(label))
  part_label  = label[idxfilled]
  part_ulabel = unique(part_label)
  part_data   = pX[idxfilled,]
  #   3. build neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  matS      = nbdstruct$mask*1.0
  matL      = (diag(rowSums(matS))-matS)
  #   4. per-class and overall : mean vectors
  meanvectors   = ammc_meanvec(part_data, part_label, part_ulabel)
  mean_Overall  = meanvectors$overall
  mean_PerClass = meanvectors$class
  #   5. adaptive scatter matrices
  adaSb = ammc_adaSb(mean_PerClass, a)
  adaSw = ammc_adaSw(part_data, part_label, b)

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION
  costS = (adaSb - (lambda*adaSw) - beta*(t(pX)%*%matL%*%pX))
  projection = aux.adjprojection(RSpectra::eigs(costS, ndim)$vectors)

  #------------------------------------------------------------------------
  ## RETURN THE RESULTS
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
