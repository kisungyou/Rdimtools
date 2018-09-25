#' Spectral Multidimensional Scaling
#'
#' \code{do.spmds} transfers the classical multidimensional scaling problem into
#' the data spectral domain using Laplace-Beltrami operator. Its flexibility
#' to use subsamples and spectral interpolation of non-reference data enables relatively
#' efficient computation for large-scale data.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param neigs number of eigenvectors to be used as \emph{spectral dimension}.
#' @param ratio percentage of subsamples as reference points.
#' @param preprocess an additional option for preprocessing the data.
#' Default is \code{"null"}. See also \code{\link{aux.preprocess}} for more details.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param symmetric one of \code{"intersect"}, \code{"union"} or \code{"asymmetric"} is supported. Default is \code{"union"}. See also \code{\link{aux.graphnbd}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Replicate the numerical example from the paper
#' #  Data Preparation
#' dim.true  = 3     # true dimension
#' dim.embed = 100   # embedding space (high-d)
#' npoints   = 1000  # number of samples to be generated
#'
#' v     = matrix(runif(dim.embed*dim.true),ncol=dim.embed)
#' coeff = matrix(runif(dim.true*npoints),  ncol=dim.true)
#' X     = coeff%*%v
#'
#' # see the effect of neighborhood size
#' out1  = do.spmds(X, neigs=100, type=c("proportion",0.1))
#' out2  = do.spmds(X, neigs=100, type=c("proportion",0.25))
#' out3  = do.spmds(X, neigs=100, type=c("proportion",0.50))
#'
#' # visualize the results
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1],out1$Y[,2],main="10% neighborhood")
#' plot(out2$Y[,1],out2$Y[,2],main="25% neighborhood")
#' plot(out3$Y[,1],out3$Y[,2],main="50% neighborhood")
#' }
#'
#' @references
#' \insertRef{aflalo_spectral_2013}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_SPMDS
#' @export
do.spmds <- function(X, ndim=2, neigs=max(2,nrow(X)/10), ratio=0.1,
                     preprocess=c("null","center","scale","cscale","decorrelate","whiten"),
                     type=c("proportion",0.1), symmetric=c("union","intersect","asymmetric")){
  #------------------------------------------------------------------------
  ## PREPROCESSING : PARAMETERS
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim and neigs
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.spmds : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  #   3. sample ratio
  if ((length(ratio)>1)||    (!((ratio<1)&&(ratio>0)))){
    stop("* do.spmds : 'ratio' is invalid.")
  }
  m_e  = as.integer(neigs)
  if ((!is.numeric(m_e))||(m_e<2)||(m_e>=n)||(length(m_e)>1)||(is.infinite(m_e))||is.na(m_e)){
    stop("* do.spmds : 'neigs' is not valid.")
  }
  #   3. preprocessing of data
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   4. neighborhood information
  nbdtype = type
  if (missing(symmetric)){
    nbdsymmetric = "union"
  } else {
    nbdsymmetric = match.arg(symmetric)
  }
  #   5. bandwidth
  # bandwidth = sum(bandwidth)
  # if (!is.numeric(bandwidth)|(bandwidth<0)|is.infinite(bandwidth)){
  #   stop("* do.spmds : 'bandwidth' should be a real number >= 0.")
  # }
  #   6. other parameters
  m_s = as.integer(n*ratio)

  #------------------------------------------------------------------------
  ## PREPROCESSING : DATA
  # 1. data preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  # 2. neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)

    # 3. construct graph laplacian : no need for scaling
  W = 1.0*nbdstruct$mask
  if (!isSymmetric(W)){
    W = (W+t(W))/2
  }
  L = diag(rowSums(W))-W
  #
  # W = exp(-(sD^2)/bandwidth) # weighted graph [0,1]
  # diag(W) = 0.0
  # if (!isSymmetric(W)){
  #   W = (W + t(W))/2
  # }
  # L = diag(rowSums(W))-W     # graph laplacian for entire dataset

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. eigendecomposition of L
  eigL = RSpectra::eigs_sym(L, m_e, which="SM") # decreasing order, top ones
  E    = eigL$vectors[,m_e:1]
  V    = diag(eigL$values[m_e:1])
  #   2. index of random samples
  indice = sample(1:n, m_s, replace=FALSE)
  #   3. pairwise squared geodesic distances for selected points
  X_ind = t(pX[indice,])
  x2    = matrix(colSums(X_ind^2),nrow=1)
  D     = array(1,c(m_s,1))%*%x2 -2*t(X_ind)%*%X_ind + t(x2)%*%array(1,c(1,m_s))

  #   4. compute
  Ylow = do.spmds.approxMDS(D, indice, E, V, ndim)

  #------------------------------------------------------------------------
  ## RETURN
  result   = list()
  result$Y = Ylow
  result$trfinfo = trfinfo
  return(result)
}


#   -----------------------------------------------------------------------
#  D      : (m.s-by-m.s) squared pairwise distances
#  sample : vector of m.s indices
#  Phi    : (M-by-K) eigenvector matrix of the set of points
#  lambda : (K-by-K) eigenvalues diagonal matrix
#  dim_o  : dimension of the embedding
#' @keywords internal
#' @noRd
do.spmds.approxMDS <- function(D, sample, Phi, lambda, dim_o){
  D = (D+t(D))/2
  sigma = 1
  n = nrow(Phi)
  N = ncol(Phi)
  Aeq = array(0,c(length(sample), n))
  for (i in 1:length(sample)){
    Aeq[i,sample[i]] = 1
  }
  Aeq_tilde = Aeq%*%Phi
  B = sigma*(t(Aeq_tilde)%*%Aeq_tilde)

  LHS = lambda[1:N,1:N]+B
  RHS = t(Aeq_tilde)
  Mat_interp = sigma*(base::solve(LHS,RHS))

  Alpha = Mat_interp%*%D%*%t(Mat_interp)

  J_Phi = Phi-(1/n)*matrix(rep(colSums(Phi),times=n),nrow=n,byrow=TRUE) ## risky
  svdJ_Phi = base::svd(J_Phi)
  S = svdJ_Phi$u
  U = diag(svdJ_Phi$d)
  V = svdJ_Phi$v
  Q1 = -U%*%t(V)%*%Alpha%*%V%*%U # coinciding with the definition ?
  eigQ1 = eigen(Q1)
  S2 = eigQ1$vectors

  eigvals = eigQ1$values # adjust for negative ones
  eigvals[(eigvals<=sqrt(.Machine$double.eps))] = sqrt(.Machine$double.eps)

  Q = (1/sqrt(2))*S%*%S2%*%diag(sqrt(eigvals))
  return(Q[,1:dim_o])
}
