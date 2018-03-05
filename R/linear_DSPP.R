#' Discriminative Sparsity Preserving Projection
#'
#' Discriminative Sparsity Preserving Projection (DSPP) is a supervised dimension reduction method
#' that employs sparse representation model to adaptively build both intrinsic adjacency graph and
#' penalty graph. It follows an integration of global within-class structure into manifold learning
#' under exploiting discriminative nature provided from label information.
#'
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param lambda regularization parameter for constructing sparsely weighted network.
#' @param rho a parameter for balancing the local and global contribution.
#'
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
#' ## generate data of 2 types with clear difference
#' diff = 15
#' dt1  = aux.gensamples(n=123)-diff;
#' dt2  = aux.gensamples(n=123)+diff;
#'
#' ## merge the data and create a label correspondingly
#' Y      = rbind(dt1,dt2)
#' label  = c(rep(1,123), rep(2,123))
#'
#' ## try different rho values
#' out1 <- do.dspp(Y, label, ndim=2, rho=0.01)
#' out2 <- do.dspp(Y, label, ndim=2, rho=0.1)
#' out3 <- do.dspp(Y, label, ndim=2, rho=1)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="rho=0.01")
#' plot(out2$Y[,1], out2$Y[,2], main="rho=0.1")
#' plot(out3$Y[,1], out3$Y[,2], main="rho=1")
#' }
#'
#' @references
#' \insertRef{gao_discriminative_2015}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_DSPP
#' @export
do.dspp <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"),
                    lambda=1.0, rho=1.0){
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
      stop("* do.dspp : no degerate class of size 1 is allowed.")
    }
  }
  if (any(is.na(label))||(any(is.infinite(label)))){stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.dspp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. lambda : positive real number
  lambda = as.double(lambda)
  if (!check_NumMM(lambda, 0, 1e+10, compact=FALSE)){stop("* do.dspp : 'lambda' is a positive real number.")}
  #   6. rho    : nonnegative real number
  rho = as.double(rho)
  if (!check_NumMM(rho, 0, 1e+10, compact=TRUE)){stop(" do.dspp : 'rho' is a nonnegative real number.")}

  #------------------------------------------------------------------------
  ## COMPUTATION 1. Preliminary Computation
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #------------------------------------------------------------------------
  ## COMPUTATION 2. Three Important Matrices
  #   1. W : adjacency, intrinsic
  W = array(0,c(n,n))
  for (i in 1:n){
    #   1-1. target vector
    xi = matrix(pX[i,],ncol=1)
    #   1-2. target matrix
    tgtidx = setdiff(which(label==label[i]),i)
    Xi = t(pX[tgtidx,])
    #   1-3. compute wi and assign it as column vector
    W[tgtidx,i] = dspp_compute_wi(xi,Xi,lambda)
  }
  #   2. B : adjacency, penalty
  #   2-1. compute radius
  radius = rep(0,n)
  for (i in 1:n){
    tgtrow  = as.vector(W[i,])
    idxmax  = which(tgtrow==max(tgtrow))
    if (length(idxmax)>1){
      idxmax = as.integer(idxmax[1])
    }
    maxdist = as.vector(pX[i,])-as.vector(pX[idxmax,])
    radius[i] = sqrt(sum(maxdist*maxdist))
  }
  #   2-2. let's iterate
  B = array(0,c(n,n))
  for (i in 1:n){
    veci = as.vector(pX[i,])
    for (j in 1:n){
      vecj = as.vector(pX[j,])
      if (i!=j){
        if (dspp_distnorm(veci,vecj)<=radius[i]){
          if (label[i]!=label[j]){
            B[i,j] = 1
          }
        }
      }
    }
  }
  #   3. compute class-wise mean
  classmean = array(0,c(length(ulabel), p))
  for (i in 1:length(ulabel)){
    partidx = which(label==ulabel[i])
    classmean[i,] = colMeans(pX[partidx,])
  }
  #   4. Sw : does not have the name on it.
  Sw = array(0,c(p,p))
  for (i in 1:n){
    #   4-1. select one vector
    cvec = pX[i,]
    clabel = label[i]
    #   4-2. select target
    tgtmat = pX[-i,]
    tgtlabel = label[-i]
    #   4-3. compute inner summation
    Smat = dspp_compute_Sw(cvec,clabel,tgtmat,tgtlabel,ulabel,classmean,as.double(radius[i]))
    Sw   = Sw+Smat
  }
  #------------------------------------------------------------------------
  ## COMPUTATION 3. set geigen setting
  #   0. symmetrize
  W = (W+t(W))/2
  B = (B+t(B))/2
  Sw= (Sw+t(Sw))/2
  #   1. Lw, Lb
  Lw = diag(rowSums(W))-W
  Lb = diag(rowSums(B))-B
  #   2. LHS and RHS
  LHS = ((t(pX)%*%Lw%*%pX) + (rho*Sw))
  RHS = (t(pX)%*%Lb%*%pX)
  #   3. use with SMALLEST ones
  projection = aux.geigen(LHS,RHS,ndim,maximal=FALSE)

  #------------------------------------------------------------------------
  ## RETURN OUTPUT

  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}



#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
dspp_compute_wi <- function(xi, Xi, lambda){
  n   = ncol(Xi)
  si  = CVXR::Variable(n)
  obj = (CVXR::p_norm(xi-(Xi%*%si),p=2) + lambda*CVXR::p_norm(si,p=1))
  constr = list(si>=0)
  prob   = CVXR::Problem(Minimize(obj), constr)
  result = solve(prob)
  return(as.vector(result$getValue(si)))
}
#' @keywords internal
#' @noRd
dspp_distnorm <- function(vec1, vec2){
  vecdiff = as.vector(vec1)-as.vector(vec2)
  output  = as.double(sqrt(sum(vecdiff*vecdiff)))
  return(output)
}
#' @keywords internal
#' @noRd
dspp_compute_Sw <- function(cvec, clabel, tgtmat, tgtlabel, ulabel, classmean, radius){
  p = length(cvec)
  Smat = array(0,c(p,p))
  ntgt = nrow(tgtmat)
  for (i in 1:ntgt){
    # 1. within radius one
    vecdiff = as.vector(cvec)-as.vector(tgtmat[i,])
    if (sqrt(sum(vecdiff*vecdiff))<radius){
      # 2. different class only
      if (tgtlabel[i]!=clabel){
        tgtvec  = (tgtmat[i,])
        tgtidx  = which(ulabel==clabel)
        meanvec = (classmean[tgtidx,])
        mvdiff  = as.vector(tgtvec)-as.vector(meanvec)
        # 3. update Smat
        Smat = Smat + outer(mvdiff,mvdiff)
      }
    }
  }
  return(Smat)
}











