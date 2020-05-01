#' (Classical) Multidimensional Scaling
#'
#' \code{do.mds} performs a classical Multidimensional Scaling (MDS) using
#' \code{Rcpp} and \code{RcppArmadillo} package to achieve faster performance than
#' \code{\link[stats]{cmdscale}}.
#'
#' @param X an \eqn{(n\times p)} matrix whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an option for preprocessing the data. Default is \code{"center"}.
#' See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#'
#' @examples
#' ## use iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' lab   = as.factor(iris[subid,5])
#'
#' ## try different preprocessing
#' out1 <- do.mds(X,ndim=2)
#' out2 <- do.mds(X,ndim=2,preprocess="cscale")
#' out3 <- do.mds(X,ndim=2,preprocess="whiten")
#'
#' ## extract embeddings for each procedure
#' Y1 <- out1$Y; Y2 <- out2$Y; Y3 <- out3$Y
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(Y1, pch=19, col=lab, main="MDS::center")
#' plot(Y2, pch=19, col=lab, main="MDS::decorrelate")
#' plot(Y3, pch=19, col=lab, main="MDS::whiten")
#' par(opar)
#'
#' @references
#' \insertRef{kruskal_multidimensional_1964}{Rdimtools}
#'
#' @concept linear_methods
#' @export
#' @rdname linear_MDS
#' @author Kisung You
do.mds <- function(X,ndim=2,preprocess=c("center","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  # Preprocessing
  if (!is.matrix(X)){stop("* do.mds : 'X' should be a matrix.")}
  myndim = round(ndim)
  myprep = ifelse(missing(preprocess), "center", match.arg(preprocess))

  #------------------------------------------------------------------------
  # Version 2 update
  output = dt_mds(X, myndim, myprep)
  return(output)

  # # 1. typecheck is always first step to perform.
  # aux.typecheck(X)
  #
  # # 2. Setting
  # #   2-1. preprocessing : center, decorrelate, or whiten
  # if (missing(preprocess)){
  #   algpreprocess = "center"
  # } else {
  #   algpreprocess = match.arg(preprocess)
  # }
  # tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  # trfinfo = tmplist$info
  # tpX     = t(tmplist$pX)
  #
  #
  # #   2-2. ndim case
  # if (!is.numeric(ndim)||(floor(ndim)>nrow(tpX))||(ceiling(ndim)<1)){
  #   stop("* do.mds : 'ndim' should have an integer value between [1,#(covariates)]")
  # }
  # tgtdim = as.integer(ndim)
  #
  # # 3. run
  # output  = method_mds(tpX);
  # eigvals = as.vector(output$eigval)
  # eigvecs = as.matrix(output$eigvec)
  #
  # # 4. result
  # result = list()
  # result$Y = t(diag(eigvals[1:tgtdim]) %*% t(eigvecs[,1:tgtdim]))
  # trfinfo$algtype   = "linear"
  # result$trfinfo    = trfinfo
  #
  # LHS = tpX %*% (tmplist$pX)
  # RHS = tpX %*% result$Y
  # result$projection = solve(LHS,RHS)
  # return(result)
}

# call for later use ------------------------------------------------------
#' @keywords internal
pydo_mds <- function(myX, mydim, myproc){
  return(do.mds(myX, ndim=mydim, preprocess=myproc))
}
