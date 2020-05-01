#' Kernel Local Discriminant Embedding
#'
#' Kernel Local Discriminant Embedding (KLDE) is a variant of Local Discriminant Embedding in that
#' it aims to preserve inter- and intra-class neighborhood information in a nonlinear manner using
#' kernel trick. \emph{Note} that the combination of kernel matrix and its eigendecomposition
#' often suffers from lacking numerical rank. For such case, our algorithm returns a warning message and
#' algorithm stops working any further due to its innate limitations of constructing weight matrix.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param t kernel bandwidth in \eqn{(0,\infty)}.
#' @param numk the number of neighboring points for k-nn graph construction.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param ktype a vector containing name of a kernel and corresponding parameters. See also \code{\link{aux.kernelcov}} for complete description of Kernel Trick.
#' @param kcentering a logical; \code{TRUE} to use centered Kernel matrix, \code{FALSE} otherwise.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#'
#' @examples
#' ## generate data of 2 types with clear difference
#' set.seed(100)
#' diff = 25
#' dt1  = aux.gensamples(n=50)-diff;
#' dt2  = aux.gensamples(n=50)+diff;
#'
#' ## merge the data and create a label correspondingly
#' X      = rbind(dt1,dt2)
#' label  = rep(1:2, each=50)
#'
#' ## try different neighborhood size
#' out1 <- do.klde(X, label, numk=5)
#' out2 <- do.klde(X, label, numk=10)
#' out3 <- do.klde(X, label, numk=20)
#'
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, pch=19, main="k=5")
#' plot(out2$Y, col=label, pch=19, main="k=10")
#' plot(out3$Y, col=label, pch=19, main="k=20")
#' par(opar)
#'
#' @references
#' \insertRef{hwann-tzongchen_local_2005}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_KLDE
#' @concept nonlinear_methods
#' @export
do.klde <- function(X, label, ndim=2, t = 1.0, numk=max(ceiling(nrow(X)/10),2),
                   preprocess=c("center","scale","cscale","decorrelate","whiten"),
                   ktype=c("gaussian",1.0), kcentering=TRUE){
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
      stop("* do.klde : no degerate class of size 1 is allowed.")
    }
  }
  if (any(is.na(label))||(any(is.infinite(label)))){stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.klde : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. t
  t = as.double(t)
  if (!check_NumMM(t,0,1e+10,compact=TRUE)){stop("* do.klde : 't' should be a positive real number.")}
  #   5. numk
  numk = as.integer(numk)
  if (!check_NumMM(numk,1,n/2,compact=FALSE)){stop("* do.klde : 'numk' should be an integer in [2,nrow(X)/2).")}
  #   6. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   Pre. preprocessing of data matrix
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   1. construct G1 (original G) and G2 (G') for same-class and different-class connectivty
  #   1-1. find k-neighborhood graph
  nbdtype   = c("knn",numk)
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric="union")
  Dmask     = nbdstruct$mask
  #   1-2. logical based-on class information
  conn_same = klde_perclass_logical(label)
  conn_diff = 1-conn_same
  #   1-3. connectivity
  G1 = Dmask*conn_same
  G2 = Dmask*conn_diff

  #   2. build AFFINITY matrix
  expD = exp(-(as.matrix(dist(pX))^2)/t)
  W1 = expD*G1
  W2 = expD*G2

  #   3. Want To Find Embedding
  #   3-1. compute kernel matrix K
  Ks = aux.kernelcov(pX,ktype)
  if (kcentering){
    K = Ks$Kcenter
  } else {
    K = Ks$K
  }

  #   3-2. LHS and RHS
  LHS = K%*%(diag(rowSums(W2))-W2)%*%K
  RHS = K%*%(diag(rowSums(W1))-W1)%*%K

  #   4. compute pseudoprojection Matrix : use top ones
  pseudoprojection = tryCatch(aux.geigen(LHS, RHS, ndim, maximal=TRUE), error=function(e)e)
  if (inherits(pseudoprojection,"error")){
    warning("* do.klde : eigendecomposition on the kernel matrix failed.")
    return(0)
  }


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = K%*%pseudoprojection
  result$trfinfo = trfinfo
  return(result)
}

#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
klde_perclass_logical <- function(label){
  n = length(label)
  out1 = matrix(0,nrow=n,ncol=n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (label[i]==label[j]){
        out1[i,j] = 1
        out1[j,i] = 1
      }
    }
  }
  diag(out1) = 1
  return(out1)
}
