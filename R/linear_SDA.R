#' Semi-Supervised Discriminant Analysis
#'
#' Semi-Supervised Discriminant Analysis (SDA) is a linear dimension reduction method
#' when label is partially missing, i.e., semi-supervised. The labeled data
#' points are used to maximize the separability between classes while
#' the unlabeled ones to estimate the intrinsic structure of the data.
#' Regularization in case of rank-deficient case is also supported via an \eqn{\ell_2}
#' scheme via \code{beta}.
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param alpha balancing parameter between model complexity and empirical loss.
#' @param beta Tikhonov regularization parameter.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
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
#' ## compare true case with missing-label case
#' out1 = do.sda(X, label)
#' out2 = do.sda(X, label_missing)
#'
#' ## visualize
#' par(mfrow=c(1,2))
#' plot(out1$Y[,1], out1$Y[,2], main="true projection")
#' plot(out2$Y[,1], out2$Y[,2], main="20% missing labels")
#'
#' @references
#' \insertRef{cai_semisupervised_2007}{Rdimtools}
#'
#' @rdname linear_SDA
#' @author Kisung You
#' @export
do.sda <- function(X, label, ndim=2, type=c("proportion",0.1), alpha=1.0, beta=1.0){
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
  if (all(!is.na(ulabel))){
    message("* Semi-Supervised Learning : there is no missing labels. Consider using Supervised methods.")
  }
  labelorder = order(label)
  labelrank  = rank(label)
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.sda : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. alpha : balancing
  alpha = as.double(alpha)
  if (!check_NumMM(alpha,0,Inf,compact=FALSE)){stop("* do.sda : 'alpha' needs to be a positive real number.")}
  #   5. beta : regularization
  beta = as.double(beta)
  if (!check_NumMM(beta,0,Inf,compact=TRUE)){stop("* do.sda : 'beta; needs to be a nonnegative real number.")}
  #   6. neighborhood type
  nbdtype = type

  #   (Implicit) preprocessing
  algpreprocess = "center"
  #   (Implicit) neighborhood symmetric
  nbdsymmetric = "union"

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing with re-labeling of data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pXoriginal   = tmplist$pX
  pX      = pXoriginal[labelorder,]
  trfinfo$algtype = "linear"
  label   = label[labelorder]
  #   2. neighborhood graph
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask
  #   3. S : binary adjacency
  S = nbdmask*1.0
  L = diag(rowSums(S))-S
  #   4. W : Weight Matrix
  idxmaxlabeled = sum(!is.na(label))
  Wl = sda_build_Wl(label[1:idxmaxlabeled])
  W  = array(0,c(n,n))
  W[1:idxmaxlabeled, 1:idxmaxlabeled] = Wl
  #   5. Itilde
  Itilde = array(0,c(n,n))
  Itilde[1:idxmaxlabeled, 1:idxmaxlabeled] = diag(idxmaxlabeled)

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. setup
  LHS = t(pX)%*%W%*%pX
  RHS = t(pX)%*%(    Itilde + (alpha*L) + (beta*diag(n)) )%*%pX
  #   2. top eigenvectors
  projection = aux.geigen(LHS, RHS, ndim, maximal=TRUE)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pXoriginal%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}




#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
sda_build_Wl <- function(vec){
  uvec = unique(vec)
  c = length(uvec)
  l = length(vec)
  output = array(0,c(l,l))
  start = 1
  for (i in 1:c){
    li = sum(vec==uvec[i])
    output[start:(start+li-1),start:(start+li-1)] = 1/li
    start = start + li
  }
  return(output)
}







