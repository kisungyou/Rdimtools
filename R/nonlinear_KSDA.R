#' Kernel Semi-Supervised Discriminant Analysis
#'
#' Kernel Semi-Supervised Discriminant Analysis (KSDA) is a nonlinear variant of
#' SDA (\code{\link{do.sda}}). For simplicity, we enabled heat/gaussian kernel only.
#' Note that this method is \emph{quite} sensitive to choices of
#' parameters, \code{alpha}, \code{beta}, and \code{t}. Especially when data
#' are well separated in the original space, it may lead to unsatisfactory results.
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
#' @param t bandwidth parameter for heat kernel.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @seealso \code{\link{do.sda}}
#' @examples
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=33)-100
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+100
#'
#' ## merge the data and create a label correspondingly
#' Y      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## copy a label and let 10% of elements be missing
#' nlabel = length(label)
#' nmissing = round(nlabel*0.10)
#' label_missing = label
#' label_missing[sample(1:nlabel, nmissing)]=NA
#'
#' ## compare true case with missing-label case
#' out1 = do.ksda(Y, label, beta=0, t=0.1)
#' out2 = do.ksda(Y, label_missing, beta=0, t=0.1)
#'
#' ## visualize
#' par(mfrow=c(1,2))
#' plot(out1$Y[,1], out1$Y[,2], main="true projection")
#' plot(out2$Y[,1], out2$Y[,2], main="20% missing labels")
#'
#' @references
#' \insertRef{cai_semisupervised_2007}{Rdimtools}
#'
#' @rdname nonlinear_KSDA
#' @author Kisung You
#' @export
do.ksda <- function(X, label, ndim=2, type=c("proportion",0.1), alpha=1.0, beta=1.0, t=1.0){
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
  if (!check_ndim(ndim,p)){
    stop("* do.ksda : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  #   4. alpha : balancing
  alpha = as.double(alpha)
  if (!check_NumMM(alpha,0,Inf,compact=FALSE)){stop("* do.ksda : 'alpha' needs to be a positive real number.")}
  #   5. beta : regularization
  beta = as.double(beta)
  if (!check_NumMM(beta,0,Inf,compact=TRUE)){stop("* do.ksda : 'beta; needs to be a nonnegative real number.")}
  #   6. neighborhood type
  nbdtype = type
  #   7. t : kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t, 0, 1e+10, compact=FALSE)){stop("* do.ksda : 't' is a bandwidth parameter for gaussian kernel.")}

  #   (Implicit) preprocessing
  algpreprocess = "center"
  #   (Implicit) neighborhood symmetric
  nbdsymmetric = "union"

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing with re-labeling of data
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX[labelorder,]
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
  #   6. computation
  K = exp(-(as.matrix(dist(pX))^2)/(2*(t^2)))

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. setup
  LHS = K%*%W%*%K
  RHS = K%*%(    Itilde + (alpha*L) + (beta*diag(n)) )%*%K
  #   2. top eigenvectors
  pseudoproj = aux.geigen(LHS, RHS, ndim, maximal=TRUE)
  #   3. find projected ones : recovering its order will be performed later.
  Y = K%*%pseudoproj

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = Y[labelrank,]
  result$trfinfo = trfinfo
  return(result)





}
