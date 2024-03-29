#' Discriminant Neighborhood Embedding
#'
#' Discriminant Neighborhood Embedding (DNE) is a supervised subspace learning method.
#' DNE tries to move multi-class data points in high-dimensional space in accordance with
#' local intra-class attraction and inter-class repulsion.
#'
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param numk the number of neighboring points for k-nn graph construction.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## load iris data
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150,50)
#' X     = as.matrix(iris[subid,1:4])
#' label = as.factor(iris[subid,5])
#'
#' ## try different numbers for neighborhood size
#' out1 = do.dne(X, label, numk=5)
#' out2 = do.dne(X, label, numk=10)
#' out3 = do.dne(X, label, numk=20)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, main="DNE::nbd size=5",  col=label, pch=19)
#' plot(out2$Y, main="DNE::nbd size=10", col=label, pch=19)
#' plot(out3$Y, main="DNE::nbd size=20", col=label, pch=19)
#' par(opar)
#'
#' @references
#' \insertRef{zhang_discriminant_2006}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_DNE
#' @concept linear_methods
#' @export
do.dne <- function(X, label, ndim=2, numk=max(ceiling(nrow(X)/10),2),
                   preprocess=c("center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.dne : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. numk
  numk = as.integer(numk)
  if (!check_NumMM(numk,1,n/2,compact=FALSE)){stop("* do.dne : 'numk' should be an integer in [2,nrow(X)/2).")}
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }
  #   5. label
  #   2. label : check and return a de-factored vector
  #   For this example, there should be no degenerate class of size 1.
  label  = check_label(label, n)
  ulabel = unique(label)
  if (any(is.na(label))||(any(is.infinite(label)))){stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #   2. neighborhood creation
  nbdtype = c("knn",numk)
  nbdsymmetric = "union"
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask
  #   3. adjacency matrix F
  matF = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      # for non-zeros, it needs to be neighbors
      if (nbdmask[i,j]==TRUE){
        if (label[i]==label[j]){
          matF[i,j] = 1.0
          matF[j,i] = 1.0
        } else {
          matF[i,j] = -1.0
          matF[j,i] = -1.0
        }
      }
    }
  }
  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. cost function :S-F
  # arbitrary regularization
  # matS = diag(rowSums(matF))-matF
  matS    = diag(rowSums(matF))
  costobj = t(pX)%*%(matS-matF)%*%pX

  #   2. minimal eigenvectors
  projection = aux.adjprojection(RSpectra::eigs(costobj, ndim, which="SR")$vectors)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
