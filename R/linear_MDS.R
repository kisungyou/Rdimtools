#' (Classical) Multidimensional Scaling
#'
#' \code{do.mds} performs a classical Multidimensional Scaling (MDS) using
#' \code{Rcpp} and \code{RcppArmadillo} package to achieve faster performance than
#' \code{\link[stats]{cmdscale}}.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an option for preprocessing the data. This supports three methods,
#' ``center'',``decorrelate'', or ``whiten''. See also \code{\link{aux.preprocess}}
#' for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#'
#'@examples
#'# generate data
#'X <- rbind(matrix(rnorm(100),nr=10),matrix(rnorm(100),nr=10)+10)
#'
#'## 1. projection onto 2 dimension.
#'output <- do.mds(X,ndim=2)
#'plot(output$Y[,1],output$Y[,2])
#'
#'## 2. different preprocessing leads to different results
#'output2 <- do.mds(X,ndim=2,preprocess="decorrelate")
#'output3 <- do.mds(X,ndim=2,preprocess="whiten")
#'par(mfrow=c(1,3))
#'plot(output$Y[,1],output$Y[,2],main="center")
#'plot(output2$Y[,1],output2$Y[,2],main="decorrelate")
#'plot(output3$Y[,1],output3$Y[,2],main="whiten")
#'
#'@references
#'\insertRef{kruskal_multidimensional_1964}{Rdimtools}
#'
#' @export
#' @rdname linear_MDS
#' @author Kisung You
do.mds <- function(X,ndim=2,preprocess="center"){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)

  # 2. Setting
  #   2-1. preprocessing : center, decorrelate, or whiten
  switch(preprocess,
         center={
           tmplist = aux.preprocess(X,type="center")
           tpX     = t(tmplist$pX)
           trfinfo = tmplist$trfinfo
         },
         decorrelate={
           tmplist = aux.preprocess(X,type="decorrelate")
           tpX     = t(tmplist$pX)
           trfinfo = tmplist$trfinfo
         },
         whiten={
           tmplist = aux.preprocess(X,type="whiten")
           tpX     = t(tmplist$pX)
           trfinfo = tmplist$trfinfo
         },
         stop("* do.mds : invalid preprocessing type. choose one of three or leave it blank.")
  )
  #   2-2. ndim case
  if (!is.numeric(ndim)||(floor(ndim)>nrow(tpX))||(ceiling(ndim)<1)){
    stop("* do.mds : 'ndim' should have an integer value between [1,#(covariates)]")
  }
  tgtdim = as.integer(ndim)

  # 3. run
  output = method_mds(tpX);
  eigvals = as.vector(output$eigval)
  eigvecs = as.matrix(output$eigvec)

  # 4. result
  result = list()
  result$Y = t(diag(eigvals[1:tgtdim]) %*% t(eigvecs[,1:tgtdim]))
  trfinfo$algtype   = "linear"
  result$trfinfo    = trfinfo

  LHS = tpX %*% (tmplist$pX)
  RHS = tpX %*% result$Y
  result$projection = solve(LHS,RHS)
  return(result)
}
