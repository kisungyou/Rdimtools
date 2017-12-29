#' Principal Component Analysis
#'
#' \code{do.pca} performs a classical principal component analysis (PCA) using
#' \code{RcppArmadillo} package for faster and efficient computation.
#'
#' A combination of \code{ndim="auto"} and \code{varratio} options is to
#' automatically decide the target dimension based on cumulative sum of
#' variance. Measured by summation of top eigenvalues from sample covariance,
#' we use the minimal summation to be larger than \code{varratio}.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension or ``auto'' option using \emph{varratio.}
#' @param cor mode of eigendecomposition. \code{FALSE} for decomposing covariance matrix,
#' and \code{TRUE} for correlation matrix.
#' @param preprocess an option for preprocessing the data. This supports three methods,
#' ``center'',``decorrelate'', or ``whiten''. See also \code{\link{aux.preprocess}}
#' for more details.
#' @param varratio a value in (0,1]. This value is only used when \code{ndim} is
#' chosen as ``auto''.
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{vars}{a vector containing variances of projected data onto principal components.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are principal components.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \dontrun{
#' # generate data
#' X <- rbind(matrix(rnorm(100),nr=10),matrix(rnorm(100),nr=10)+10)
#'
#' ## 1. projection using 2 principal components
#' output <- do.pca(X,ndim=2)
#' plot(output$Y[,1],output$Y[,2])
#'
#' ## 2. automatic detection of target dimension accounting for 98% of variance
#' output <- do.pca(X,ndim="auto",varratio=0.98)           # perform PCA
#' plot(seq_len(length(output$vars)),output$vars,type="b") # plot variances
#' }
#'
#' @author Kisung You
#' @references
#' \insertRef{pearson_liii._1901}{Rdimtools}
#'
#' @rdname linear_PCA
#' @export
do.pca <- function(X,ndim="auto",cor=FALSE,preprocess="center",varratio=0.9){
  # 1. typecheck is always first step to perform.
  aux.typecheck(X)

  # 2. preprocessing
  #   2-1. center,decorrelate, or whiten
  switch(preprocess,
         center={
           tmplist = aux.preprocess(X,"center")
           pX      = tmplist$pX
           trfinfo = tmplist$info
         },
         decorrelate={
           tmplist = aux.preprocess(X,"decorrelate")
           pX      = tmplist$pX
           trfinfo = tmplist$info
         },
         whiten={
           tmplist = aux.preprocess(X,"whiten")
           pX      = tmplist$pX
           trfinfo = tmplist$info
         },
         stop("* do.pca : invalid preprocessing type. choose one of three or leave it blank.")
  )
  #   2-2. correlation or covariance matrix
  if (cor){
    psdX = cor(pX)
  } else if (cor==FALSE){
    psdX = cov(pX)
  } else{
    stop("* do.pca : invalid 'cor' type. should be logical.")
  }

  # 3. run method_pca
  output = method_pca(psdX)
  eigvals = rev(output$eigval)
  eigvecs = output$eigvec[,rev(seq_len(length(eigvals)))]

  # 4. branching
  if (ndim=="auto"){
    if (!is.numeric(varratio)||(varratio>1)||(varratio<=0)){
      stop("* do.pca : 'varratio' should be in (0,1]")
    }
    tgtidx = as.integer(min(which((cumsum(eigvals)/sum(eigvals))>=varratio)))
    if (tgtidx<1){
      tgtidx = as.integer(1);
    }
    if (tgtidx>length(eigvals)){
      tgtidx = as.integer(length(eigvals))
    }
  } else{
    if (!is.numeric(ndim)||(ndim<1)||(ndim>length(eigvals))||is.na(ndim)){
      stop("* do.pca : 'ndim' should be in [1,#(covariates)]")
    }
    tgtidx = as.integer(ndim)
  }

  # 5. result
  partials = aux.adjprojection(eigvecs[,1:tgtidx])
  result = list()
  result$Y          = (pX %*% partials)
  result$vars       = eigvals[1:tgtidx]

  result$projection = partials
  trfinfo$algtype   = "linear"
  result$trfinfo    = trfinfo
  return(result)
}
