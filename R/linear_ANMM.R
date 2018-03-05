#' Average Neighborhood Margin Maximization
#'
#' Average Neighborhood Margin Maximization (ANMM) is a supervised method
#' for feature extraction. It aims to find a projection mapping in the following manner;
#' for each data point, the algorithm tries to pull the neighboring points in the
#' same class while pushing neighboring points of different classes far away. It is known
#' that ANMM does suffer less from small sample size problem, which is bottleneck for LDA.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param No neighborhood size for same-class data points; either a constant number or
#' a vector of length-\eqn{n} can be provided, as long as the values reside in \eqn{[2,n]}.
#' @param Ne neighborhood size for different-class data points; either a constant number or
#' a vector of length-\eqn{n} can be provided, as long as the values reside in \eqn{[2,n]}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate 3 different groups of data X and label vector
#' x1 = matrix(rnorm(4*10), nrow=10)-20
#' x2 = matrix(rnorm(4*10), nrow=10)
#' x3 = matrix(rnorm(4*10), nrow=10)+20
#' X  = rbind(x1, x2, x3)
#' label = c(rep(1,10), rep(2,10), rep(3,10))
#'
#' ## perform ANMM on different choices of neighborhood size
#' out1 = do.anmm(X, label, No=6, Ne=6)
#' out2 = do.anmm(X, label, No=2, Ne=10)
#' out3 = do.anmm(X, label, No=10,Ne=2)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="(No,Ne)=(6,6)")
#' plot(out2$Y[,1], out2$Y[,2], main="(No,Ne)=(2,10)")
#' plot(out3$Y[,1], out3$Y[,2], main="(No,Ne)=(10,2)")
#' }
#'
#'
#' @references
#' \insertRef{wang_feature_2007}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_ANMM
#' @export
do.anmm <- function(X, label, ndim=2, preprocess=c("null","center","scale","cscale","decorrelate","whiten"),
                    No=ceiling(nrow(X)/10), Ne=ceiling(nrow(X)/10)){
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
      stop("* do.anmm : no degerate class of size 1 is allowed.")
    }
  }
  if (any(is.na(label))||(any(is.infinite(label)))){
    stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")
  }

  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.anmm : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  algpreprocess = match.arg(preprocess)
  #   5. No : size for Homogeneous   Neighborhood
  #      Ne :          Heterogenoeus Neighborhood
  vecNo = anmm_nbdstructure(No, n, 1)
  vecNe = anmm_nbdstructure(Ne, n, 2)

  #------------------------------------------------------------------------
  ## COMPUTATION
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. extract 2 types of neighborhood information
  D = as.matrix(dist(pX, method="euclidean"))
  listNo = anmm_find_No(D, label, vecNo)
  listNe = anmm_find_Ne(D, label, vecNe)
  #   3. compute scatter and compactness
  S = anmm_computeSC(pX, listNe)
  C = anmm_computeSC(pX, listNo)
  #   4. extract eigenvector
  eigSC      = RSpectra::eigs_sym(S-C, ndim, which="LA")
  projection = matrix(eigSC$vectors, nrow=p)
  #   5. adjust eigenvectors
  projection = aux.adjprojection(projection)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}

# auxiliary for ANMM ------------------------------------------------------
# 1. nbdstructure input argument
#' @keywords internal
#' @noRd
anmm_nbdstructure <- function(sizevec, n, Ntype){
  tmp = as.vector(round(sizevec))
  if (length(tmp)==1){
    nbdvec = rep(tmp, n)
  } else if (length(tmp)==n){
    nbdvec = tmp
  } else {
    if (Ntype==1){
      stop("* do.anmm : homogeneous neighborhood input is invalid.")
    } else {
      stop("* do.anmm : heterogeneous neighborhood input is invalid.")
    }
  }
 if ((any(nbdvec<1))||(any(nbdvec>n))){
   if (Ntype==1){
     stop("* do.anmm : range of values from No is invalid.")
   } else {
     stop("* do.anmm : range of values from Ne is invalid.")
   }
 }
  return(nbdvec)
}

# HERE comes the real problem. Say we needed to choose 3-Ho neighborhood.
# However, it is still possible that we only have 1 element in each class.
# For this case, I already cleared it out in the preprocessin step,
# by not allowing the degenerate class of size 1. In both cases, it returns
# a list containing membership structure at each node.
#
# 2. find homogeneous neighborhood
#' @keywords internal
#' @noRd
anmm_find_No <- function(matD, label, vecNo){
  n = length(label)
  if (nrow(matD)!=n){stop("* do.anmm : I don't know why it stopped 1.")}
  output = list()
  numNo  = rep(0,n)
  for (i in 1:n){
    # compute possible values by taking minimization
    clabel  = which(label==label[i])
    nclabel = sum(label==label[i])
    nselect = round(min(nclabel, vecNo[i]))
    numNo[i]= nselect
    # find the distance
    tgtdist  = matD[i,clabel]
    smindex  = which(order(tgtdist)<=(nselect+1))
    tgtlabel = setdiff(clabel[smindex], round(i))
    output[[i]] = tgtlabel
  }

  return(output)
}

# 3. find heterogeneous neighborhood
#' @keywords internal
#' @noRd
anmm_find_Ne <- function(matD, label, vecNe){
  n = length(label)
  if (nrow(matD)!=n){stop("* do.anmm : I don't know why it stopped 2.")}
  output = list()
  numNe  = rep(0,n)
  for (i in 1:n){
    clabel  = which(label!=label[i])
    nclabel = length(clabel)
    nselect = round(min(nclabel, vecNe[i]))
    numNe[i]= nselect

    tgtdist = matD[i,clabel]
    smindex = which(order(tgtdist)<=(nselect+1))
    tgtlabel= setdiff(clabel[smindex], round(i))
    output[[i]] = tgtlabel
  }

  return(output)
}

# 4. compute Scatterness or Compactness
#' @keywords internal
#' @noRd
anmm_computeSC <- function(X, memlist){
  n = nrow(X)
  p = ncol(X)
  S = array(0, c(p,p))
  for (i in 1:n){
    xi = as.vector(X[i,])
    tgtvecs = memlist[[i]]
    tgtsize = length(tgtvecs)
    Stmp    = array(0,c(p,p))
    for (k in 1:tgtsize){
      xk   = as.vector(X[tgtvecs[k],])
      xdiff= xi-xk
      Stmp = Stmp + outer(xdiff,xdiff)
    }
    Stmp = Stmp/tgtsize
    S    = S + Stmp
  }
  return(S)
}
