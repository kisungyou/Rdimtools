#' Robust Principal Component Analysis
#'
#' Robust PCA (RPCA) is not like other methods in this package as finding explicit low-dimensional embedding with reduced number of columns.
#' Rather, it is more of a decomposition method of data matrix \eqn{X}, possibly noisy, into low-rank and sparse matrices by
#' solving the following,
#' \deqn{\textrm{minimize}\quad \|L\|_* + \lambda \|S\|_1}
#' \deqn{\textrm{subject to}\quad L+S=X}
#' where \eqn{L} is a low-rank matrix, \eqn{S} is a sparse matrix and \eqn{\|\cdot\|_*} denotes nuclear norm, i.e., sum of singular values. Therefore,
#' it should be considered as \emph{preprocessing} procedure of denoising. Note that after RPCA is applied, \eqn{L} should be used
#' as kind of a new data matrix for any manifold learning scheme to be applied. \pkg{CVXR} was used for \code{do.rpca}.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param lambda parameter for the sparsity term \eqn{\|S\|_1}. Default value is given accordingly to the referred paper.
#' @param preprocess an additional option for preprocessing the data. Default is ``center'', and other methods of ``decorrelate'' and ``whiten'' are supported. See also \code{\link{aux.preprocess}} for more details.
#'
#'
#' @return a named list containing
#' \describe{
#' \item{L}{an \eqn{(n\times p)} low-rank matrix.}
#' \item{S}{an \eqn{(n\times p)} sparse matrix.}
#' \item{cvxr.status}{``optimal'' denotes the problem was solved. See\code{\link[CVXR]{psolve}} for more details on solvability.}
#' \item{cvxr.niters}{the number of iterations taken.}
#' \item{cvxr.solver}{type of solver used by \pkg{CVXR}.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Load Iris data and put some noise
#' data(iris)
#' noise = 0.2
#' X = as.matrix(iris[,1:4])
#' X = X + matrix(noise*rnorm(length(X)), nrow=nrow(X))
#'
#' ## Compare 2 methods : {PCA} vs {RPCA + PCA}
#' rpca     = do.rpca(X)
#' out_pca  = do.pca(X,      ndim=2, preprocess="center")
#' out_rpca = do.pca(rpca$L, ndim=2, preprocess="center")
#'
#' ## Visualize
#' par(mfrow=c(1,2))
#' plot(out_pca$Y[,1], out_pca$Y[,2], main="PCA")
#' plot(out_rpca$Y[,1], out_rpca$Y[,2], main="RPCA")
#' }
#'
#' @references
#' \insertRef{candes_robust_2011}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_RPCA
#' @export
do.rpca <- function(X, lambda=sqrt(1/(max(dim(X)))), preprocess=c("center","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. lambda
  if (length(as.vector(lambda))!=1){
    stop("* do.rpca : input 'lambda' should be a constant.")
  }
  lambda = as.double(lambda)
  if ((is.na(lambda))||(is.infinite(lambda))||(lambda<=0)){
    stop("* do.rpca : 'lambda' should be a POSITIVE real number greater than 0.")
  }

  #   3. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "nonlinear"

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION USING CVXR
  #   1. define variables
  L = CVXR::Variable(n,p)
  M = CVXR::Variable(n,p)
  #   2. objective and constraint
  objfunc = CVXR::Minimize(cvxr_norm(L,"nuc")+lambda*cvxr_norm(M,1))
  constr  = list(L+M==pX)
  #   3. problem define and solve
  probRPCA = CVXR::Problem(objfunc, constr)
  solvRPCA = solve(probRPCA)

  #------------------------------------------------------------------------
  ## REPORT RESULTS
  result = list()
  result$L = solvRPCA$getValue(L)
  result$S = solvRPCA$getValue(M)
  result$cvxr.status = solvRPCA$status
  result$cvxr.niters = solvRPCA$num_iters
  result$cvxr.solver = solvRPCA$solver
  return(result)
}
