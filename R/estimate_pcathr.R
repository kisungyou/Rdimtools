#' PCA Thresholding with Accumulated Variance
#'
#' Principal Component Analysis exploits sample covariance matrix whose
#' eigenvectors and eigenvalues are principal components and projected
#' variance, correspondingly. Given \code{varratio}, it thresholds the
#' accumulated variance and selects the estimated dimension. Note that other than
#' linear submanifold case, the naive selection scheme from this algorithm
#' lacks flexibility in discovering intrinsic dimension.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param varratio target explainability for accumulated variance in \eqn{(0,1)}.
#'
#' @return a named list containing containing \describe{
#' \item{estdim}{estimated dimension according to \code{varratio}.}
#' \item{values}{eigenvalues of sample covariance matrix.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate 3-dimensional normal data
#' X = matrix(rnorm(100*3), nrow=100)
#'
#' ## replicate 3 times with translations
#' Y = cbind(X-10,X,X+10)
#'
#' ## use PCA thresholding estimation with 95% variance explainability
#' ## desired return is for dimension 3.
#' output   = est.pcathr(Y)
#' pmessage = paste("* estimated dimension is ",output$estdim, sep="")
#' print(pmessage)
#'
#' ## use screeplot
#' opar <- par(no.readonly=TRUE)
#' plot(output$values, main="scree plot")
#' par(opar)
#' }
#'
#' @seealso \code{\link{do.pca}}
#' @rdname estimate_pcathr
#' @author Kisung You
#' @export
est.pcathr <- function(X, varratio=0.95){
  # ------------------------------------------------
  # PREPROCESSING
  aux.typecheck(X)
  X = as.matrix(X)
  varratio = as.double(varratio)
  if(!check_NumMM(varratio,0,1,compact=FALSE)){
    stop("* est.pcathr : 'varratio' should be a real number in (0,1).")
  }

  # ------------------------------------------------
  # MAIN COMPUTATION
  values = eigen(cov(X))$values
  accval = cumsum(values)
  maxval = max(accval)

  estdim = min(which(accval >= (varratio*maxval)))

  # ------------------------------------------------
  # RETURN RESULTS
  output = list()
  output$estdim = estdim
  output$values = values
  return(output)
}
