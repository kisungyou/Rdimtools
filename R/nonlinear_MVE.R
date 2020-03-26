#' Minimum Volume Embedding
#'
#' Minimum Volume Embedding (MVE) is a nonlinear dimension reduction
#' algorithm that exploits semidefinite programming (SDP), like MVU/SDE.
#' Whereas MVU aims at stretching through all direction by maximizing
#' \eqn{\sum \lambda_i}, MVE only opts for unrolling the top eigenspectrum
#' and chooses to shrink left-over spectral dimension. For ease of use,
#' unlike kernel PCA, we only made use of Gaussian kernel for MVE. Note that
#' we adopted \href{https://CRAN.R-project.org/package=Rcsdp}{Rcsdp} package in that
#' when given large-scale dataset, it may result in extremely deteriorated computational performance.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param knn size of \eqn{k}-nn neighborhood.
#' @param kwidth bandwidth for Gaussian kernel.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null". See also \code{\link{aux.preprocess}} for more details.
#' @param tol stopping criterion for incremental change.
#' @param maxiter maximum number of iterations allowed.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \donttest{
#' ## generate ribbon-shaped data
#' ## in order to pass CRAN pretest, n is set to be small.
#' X = aux.gensamples(dname="ribbon",n=50)
#'
#' ## Compare MVU and MVE
#' #  Note that MVE actually requires much larger number of iterations
#' #  Here, due to CRAN limit, it was set as 7.
#' outMVU5  <- do.mvu(X, ndim=2, type=c("knn",5), projtype="kpca")
#' outMVE5  <- do.mve(X, ndim=2, knn=5, maxiter=7)
#'
#' ## Visualize two comparisons
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(outMVU5$Y,  main="MVU (k=5)")
#' plot(outMVE5$Y,  main="MVE (k=5)")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{shaw_minimum_2007}{Rdimtools}
#'
#' @seealso \code{\link{do.mvu}}
#' @author Kisung You
#' @rdname nonlinear_MVE
#' @export
do.mve <- function(X, ndim=2, knn=ceiling(nrow(X)/10), kwidth=1.0,
                   preprocess=c("null","center","scale","cscale","whiten","decorrelate"),
                   tol=1e-4, maxiter=1000){
  #------------------------------------------------------------------------
  ## PARAMETER CHECK
  #   1. X : data matrix
  aux.typecheck(X)
  #   2. ndim : target dimension
  if (!check_ndim(ndim,ncol(X))){
    stop("* do.mve : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  #   3. knn    : the number for k-nn graph
  if (length(as.vector(knn))!=1){
    stop("* do.mve : knn should be a constant integer number.")
  }
  knn = as.integer(knn)
  if ((knn<1)||(knn>=nrow(X))||(is.na(knn))||(is.infinite(knn))){
    stop("* do.mve : knn should be [1,#(covariates)).")
  }
  #   4. kwidth : kernel bandwidth
  if (length(as.vector(kwidth))!=1){
    stop("* do.mve : kernel bandwidth should be a constant number.")
  }
  if ((kwidth < 0)||(is.na(kwidth))||(is.infinite(kwidth))){
    stop("* do.mve : kernel bandwidth should be a nonnegative real number.")
  }
  ktype = c("gaussian",as.double(kwidth))
  #   5. tol    : tolerance level for stopping criterion
  if (length(as.vector(tol))!=1){
    stop("* do.mve : 'tol' should be a number.")
  }
  if ((tol<=0)||(tol>=1)||(is.na(tol))||(is.infinite(tol))){
    stop("* do.mve : 'tol' should be a small positive number, possibly in (0,1).")
  }
  #   6. preprocess
  algpreprocess = match.arg(preprocess)
  #   7. maxiter
  maxiter = round(maxiter)
  if (length(as.vector(maxiter))!=1){
    stop("* do.mve : 'maxiter' should be a constant integer.")
  }
  if ((maxiter<3)||(is.na(maxiter))||(is.infinite(maxiter))){
    stop("* do.mve : 'maxiter' should be a relatively large positive integer.")
  }

  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. datapreprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="nonlinear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #   2. form affinity matrix A
  KernelMap = aux.kernelcov(pX,ktype)
  A         = KernelMap$K
  #   3. use A to find a binary connectivity matrix C via k-nearest neighbors
  nbdtype   = c("knn",knn)
  nbdstruct = aux.graphnbdD(A,type=nbdtype,symmetric="union")
  C         = nbdstruct$mask # TRUE for connected, FALSE not connected.

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  # 1. initialize K=A
  N    = nrow(A)
  Kold = A
  cvtgap = 1000
  itercount = 1
  while (cvtgap > tol){
    # 2. solve for the eigenvectors of K
    eigKold = eigen(Kold)$vectors
    B1 = eigKold[,1:ndim]
    B2 = eigKold[,(ndim+1):ncol(eigKold)]
    B  = -(B1%*%t(B1))+(B2%*%t(B2))
    # 3. solve SDP via Rcsdp
    Knew = mve_single_csdp(A,B,C)
    # 4. update the cvtgap : I will use Frobenius norm
    cvtgap = base::norm(Kold-Knew,type="F")
    # 5. update Kold -> we will use Kold forever
    Kold   = Knew
    itercount = itercount+1
    if (itercount > maxiter){
      cvtgap = tol/100;
    }
  }
  # 7. now, Kold is our kernel matrix, so apply Kernel PCA
  tY     = aux.kernelprojection(Kold, ndim)

  #------------------------------------------------------------------------
  ## RETURN OUTPUTKnew = as.matrix(solprob$getValue(Ktmp), nrow=N)
  result = list()
  result$Y = t(tY)
  result$trfinfo = trfinfo
  return(result)
}


#' @keywords internal
#' @noRd
mve_single_cvxr <- function(A, B, C){
  N = nrow(B)
  Ktmp = CVXR::Variable(N,N,PSD=TRUE)
  obj  = Maximize(matrix_trace(Ktmp%*%B))
  constr1 = list(CVXR::sum_entries(Ktmp)==0)
  constr2 = list()
  iter = 1
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (C[i,j]){  # if two nodes are connected
        constr2[[iter]] = ((Ktmp[i,i]+Ktmp[j,j]-Ktmp[i,j]-Ktmp[j,i])==(A[i,i]+A[j,j]-A[i,j]-A[j,i]))
      }
      iter = iter+1 # update iteration counter
    }
  }
  prob = CVXR::Problem(obj, c(constr1, constr2))
  solprob = solve(prob)
  Knew = as.matrix(solprob$getValue(Ktmp), nrow=N)
  return(Knew)
}

#' @keywords internal
#' @noRd
mve_single_csdp <- function(A, B, C){
  # 1. settings
  N = nrow(B)
  nconstraints = sum(C)/2
  setC = list(B)
  setK = list(type="s", size=N)
  # 2. iterating for conditions
  #   2-1. setup
  setA = list()
  setb = c()
  #   2-2. iterative conditions
  iter = 1
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (C[i,j]==TRUE){
        tmpA = list(simple_triplet_sym_matrix(i=c(i,i,j,j),j=c(i,j,i,j),v=c(1,-1,-1,1),n=N))
        tmpb = A[i,i]+A[j,j]-A[i,j]-A[j,i]
        setA[[iter]] = tmpA
        setb[iter]   = tmpb
        iter = iter+1
      }
    }
  }
  #   2-3. sum to zero
  setA[[iter]] = list(matrix(1,N,N)/N)
  setb[iter]   = 0
  outCSDP   = (csdp(setC,setA,setb,setK, csdp.control(printlevel=0)))
  return(matrix(outCSDP$X[[1]], nrow=N))
}
