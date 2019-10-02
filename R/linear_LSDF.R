#' Locality Sensitive Discriminant Feature
#'
#' Locality Sensitive Discriminant Feature (LSDF) is a semi-supervised feature selection method.
#' It utilizes both labeled and unlabeled data points in that labeled points are used to maximize
#' the margin between data opints from different classes, while labeled ones are used to discover
#' the geometrical structure of the data space.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels. It should contain \code{NA} elements for missing label.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param gamma within-class weight parameter for same-class data.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=33)-100
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+100
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
#' ## try different neighborhood sizes
#' out1 = do.lsdf(X, label_missing, type=c("proportion",0.01))
#' out2 = do.lsdf(X, label_missing, type=c("proportion",0.1))
#' out3 = do.lsdf(X, label_missing, type=c("proportion",0.25))
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="1% connectivity")
#' plot(out2$Y[,1], out2$Y[,2], main="10% connectivity")
#' plot(out3$Y[,1], out3$Y[,2], main="25% connectivity")
#' }
#'
#' @references
#' \insertRef{cai_locality_2007}{Rdimtools}
#'
#' @rdname linear_LSDF
#' @author Kisung You
#' @export
do.lsdf <- function(X, label, ndim=2, type=c("proportion",0.1),
                    preprocess=c("null","center","scale","cscale","whiten","decorrelate"), gamma=100){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label : check and return a de-factored vector
  #   For this example, there should be no degenerate class of size 1.
  if (missing(label)){
    stop("* Semi-Supervised Learning : 'label' is required. For it not provided, consider using Unsupervised methods.")
  }
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
  if (!check_ndim(ndim,p)){
    stop("* do.lsdf : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. type
  nbdtype = type
  nbdsymmetric = "union"
  #   5. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   6. gamma
  gamma = as.double(gamma)
  if (!check_NumMM(gamma,1,1e+10)){stop("* do.lsdf : 'gamma' is a large positive real number.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. build neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION FOR LSDF
  #   1. build Within- and between-class weights
  Sb = array(0,c(n,n))
  Sw = array(0,c(n,n))
  for (i in 1:(n-1)){
    class1 = label[i]
    for (j in (i+1):n){
      class2 = label[j]
      if (((!is.na(class1))&&(!is.na(class2)))&&(class1==class2)){
        Sw[i,j] = gamma
        Sw[j,i] = gamma
      } else if ((isTRUE(nbdmask[i,j])||isTRUE(nbdmask[j,i]))&&(is.na(class1)||is.na(class2))){
        Sw[i,j] = 1.0
        Sw[j,i] = 1.0
      }
      if (((!is.na(class1))&&(!is.na(class2)))&&(class1!=class2)){
        Sb[i,j] = 1.0
        Sb[j,i] = 1.0
      }
    }
  }
  #   2. laplacian graphs
  Lw = diag(rowSums(Sw))-Sw
  Lb = diag(rowSums(Sb))-Sb
  #   3. compute feature scores
  fscore = rep(0,p)
  for (j in 1:p){
    fr = as.vector(pX[,j])
    term1 = sum(as.vector(Lb%*%matrix(fr))*fr)
    term2 = sum(as.vector(Lw%*%matrix(fr))*fr)
    fscore[j] = term1/term2
  }
  #   4. find the largest ones
  idxvec = base::order(fscore, decreasing=TRUE)[1:ndim]
  #   5. find the projection matrix
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$featidx = idxvec
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
