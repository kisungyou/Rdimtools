#' Kernel Marginal Fisher Analysis
#'
#' Kernel Marginal Fisher Analysis (KMFA) is a nonlinear variant of MFA using kernel tricks.
#' For simplicity, we only enabled a heat kernel of a form
#' \deqn{k(x_i,x_j)=\exp(-d(x_i,x_j)^2/2*t^2)}
#' where \eqn{t} is a bandwidth parameter. Note that the method is far sensitive to the choice of \eqn{t}.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param k1 the number of same-class neighboring points (homogeneous neighbors).
#' @param k2 the number of different-class neighboring points (heterogeneous neighbors).
#' @param t bandwidth parameter for heat kernel in \eqn{(0,\infty)}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
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
#' ## try different numbers for neighborhood size
#' out1 = do.kmfa(X, label, k1=5, k2=5, t=1)
#' out2 = do.kmfa(X, label, k1=5, k2=5, t=2)
#'
#' ## visualize
#' opar = par(mfrow=c(1,2), no.readonly=TRUE)
#' plot(out1$Y, main="bandwidth=1")
#' plot(out2$Y, main="bandwidth=2")
#' par(opar)
#' }
#' @references
#' \insertRef{yan_graph_2007}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_KMFA
#' @export
do.kmfa <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"),
                   k1=max(ceiling(nrow(X)/10),2), k2=max(ceiling(nrow(X)/10),2), t=1.0){
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
  if (any(is.na(label))||(any(is.infinite(label)))){warning("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.kmfa : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }
  #   5. k1 and k2
  k1 = as.integer(k1)
  k2 = as.integer(k2)
  if (!check_NumMM(k1,1,n/2,compact=FALSE)){stop("* do.kmfa : 'k1' should be an integer in [2,nrow(X)/2).")}
  if (!check_NumMM(k2,1,n/2,compact=FALSE)){stop("* do.kmfa : 'k2' should be an integer in [2,nrow(X)/2).")}
  #   6. t : kernel bandwidth
  t = as.double(t)
  if (!check_NumMM(t, 0, 1e+10, compact=FALSE)){stop("* do.kmfa : 't' is a bandwidth parameter for gaussian kernel.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #   2. PCA preprocessing
  eigtest = eigen(cov(pX), only.values=TRUE)
  pcadim  = sum(eigtest$values > 0)
  if (pcadim <= ndim){
    warning("* do.lspp : target 'ndim' is larger than intrinsic data dimension achieved from PCA.")
    projection_first = diag(p)
    pcapX = pX
  } else{
    projection_first = aux.adjprojection(eigen(cov(pX))$vectors[,1:pcadim])
    pcapX = pX%*%projection_first
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR KMFA
  #   1. we need kernel matrix K first
  K = exp(-(as.matrix(dist(pcapX))^2)/(2*(t^2)))
  #   2. nbdlogical part
  #   2-1. nbdlogical is modified
  logicalmat = nbdlogical_kmf(K, label, k1, k2)
  mat_hom    = logicalmat$hom
  mat_het    = logicalmat$het
  #   2-2. OR class under symmetrization
  W  = array(as.logical(mat_hom + t(mat_hom)),c(n,n))*1.0; diag(W)=0
  Wp = array(as.logical(mat_het + t(mat_het)),c(n,n))*1.0; diag(Wp)=0
  #   2-3. D and Dp as well
  Dp = diag(rowSums(Wp))
  D  = diag(rowSums(W))
  #   3. formulation as generalized eigenvalue problem
  LHS = K%*%(D-W)%*%K
  RHS = K%*%(Dp-Wp)%*%K
  #   4. solve for alphas : use lowest
  alphas = aux.geigen(LHS, RHS, ndim, maximal=FALSE)

  #------------------------------------------------------------------------
  ## RETURN THE RESULTS
  #   1. compute naive projection
  projected = K%*%alphas
  #   2. modify by column
  for (i in 1:ndim){
    alpha = as.vector(alphas[,i])
    Kalpha = as.vector(K%*%as.matrix(alpha))
    lambda = sqrt(alpha*Kalpha)
    projected[,i] = projected[,i]/lambda
  }
  #   3. return adjusted
  result = list()
  result$Y = projected
  result$trfinfo = trfinfo
  return(result)
}







#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
nbdlogical_kmf <- function(K, label, khomo, khet){
  n = nrow(K)
  D = array(0,c(n,n))
  for (i in 1:n){
    for (j in 1:n){
      D[i,j] = sqrt(K[i,i]+K[j,j]-2*K[i,j])
    }
  }
  # 1. homogeneous logical matrix
  logical_hom = array(FALSE,c(n,n))
  for (i in 1:n){
    # 1-1. index of same class
    idxhom = setdiff(which(label==label[i]), i)
    # 1-2. which is the smallest numk ?
    partD = as.vector(D[i,idxhom])
    # 1-3. partially smallest ones
    partidx = which(      partD <= max(sort(partD)[1:max(min(khomo, length(idxhom)),1)])    )
    # 1-4. adjust idxhom
    idxhomadj = idxhom[partidx]
    logical_hom[i,idxhomadj] = TRUE
  }
  # 2. heterogeneous logical matrix
  logical_het = array(FALSE,c(n,n))
  for (i in 1:n){
    idxhet = which(label!=label[i])
    partD = as.vector(D[i, idxhet])
    partidx = which(partD <= max(sort(partD)[1:max(min(khet, length(idxhet)),1)]))
    idxhetadj = idxhet[partidx]
    logical_het[i,idxhetadj] = TRUE
  }
  output = list()
  output$hom = logical_hom
  output$het = logical_het
  return(output)
}

