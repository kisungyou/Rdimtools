#' Constrained Graph Embedding
#'
#' Constrained Graph Embedding (CGE) is a semi-supervised embedding method that incorporates
#' partially available label information into the graph structure that find embeddings
#' consistent with the labels.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' @param label a length-\eqn{n} vector of data class labels. It should contain \code{NA} elements for missing label.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is \code{"null"}. See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,2:4])
#' label = as.integer(iris[,5])
#' lcols  = as.factor(label)
#'
#' ## copy a label and let 10% of elements be missing
#' nlabel = length(label)
#' nmissing = round(nlabel*0.10)
#' label_missing = label
#' label_missing[sample(1:nlabel, nmissing)]=NA
#'
#' ## try different neighborhood sizes
#' out1 = do.cge(X, label_missing, type=c("proportion",0.01))
#' out2 = do.cge(X, label_missing, type=c("proportion",0.1))
#' out3 = do.cge(X, label_missing, type=c("proportion",0.25))
#'
#' ## visualize
#' opar = par(mfrow=c(1,3), pty="s")
#' plot(out1$Y[,1], out1$Y[,2], main="1% connected", pch=19, col=lcols)
#' plot(out2$Y[,1], out2$Y[,2], main="10% connected", pch=19, col=lcols)
#' plot(out3$Y[,1], out3$Y[,2], main="25% connected", pch=19, col=lcols)
#' par(opar)
#'
#' @references
#' \insertRef{he_graph_2009}{Rdimtools}
#'
#' @rdname nonlinear_CGE
#' @author Kisung You
#' @export
do.cge <- function(X, label, ndim=2, type=c("proportion",0.1),
                   preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label : check and return a de-factored vector
  #   For this example, there should be no degenerate class of size 1.
  if (missing(label)){
    stop("* Semi-Supervised Learning : 'label' is required. For it not provided, consider using Unsupervised methods.")
  }
  label  = check_label(label, n)
  ulabel = unique(label)
  if (all(!is.na(ulabel))){
    message("* Semi-Supervised Learning : there is no missing labels. Consider using Supervised methods.")
  }
  if (any(is.infinite(ulabel))){
    stop("* Semi-Supervised Learning : Inf is not allowed in label.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.cge : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. type
  nbdtype = type
  nbdsymmetric = "union"
  #   5. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. re-arrange with label information
  id.labeled = which(!is.na(label))
  id.notlabd = which(is.na(label))
  p = length(label[id.labeled])

  pXnew = rbind(pX[id.labeled,], pX[id.notlabd,])  # first row section is for labeled data
  label = c(label[id.labeled], label[id.notlabd])  #
  ulabe = unique(label[1:p])
  n = ncol(pXnew)
  m = nrow(pXnew)
  c = length(ulabe)

  #   3. build neighborhood information among labeled
  nbdstruct = aux.graphnbd(pXnew,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)


  #   4. construct U matrix
  U1 = array(0,c(m,c))
  for (i in 1:c){
    idU1 = which(label==ulabel[i])
    U1[idU1,i] = rep(1,length(idU1))
  }
  U2top = array(0,c(p,m-p))
  U2bot = diag(m-p)
  U     = cbind(U1,rbind(U2top, U2bot))

  #   5. construct L and D
  W   = nbdstruct$mask*1   # should be (m x m) matrix
  D   = diag(base::rowSums(W))
  L   = D-W


  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART
  #   1. generalized eigenvalue problem; numerical error serious..
  llterm = t(U)%*%L%*%U; alpha1 = cge.minimaladd(llterm)$alpha
  rrterm = t(U)%*%D%*%U; alpha2 = cge.minimaladd(rrterm)$alpha
  glterm = llterm + max(alpha1, alpha2)*diag(nrow(llterm))
  grterm = rrterm + max(alpha1, alpha2)*diag(nrow(rrterm))

  geigs  = geigen::geigen(glterm, grterm) # increasing order
  idmin  = max(which.min(geigs$values > 10*.Machine$double.eps), 2)
  Z      = geigs$vectors[,idmin:(idmin+ndim-1)]

  #   2. reconstruct
  pY     = U%*%Z

  ########################################################################
  ## 5. return output
  result = list()
  result$Y = pY
  trfinfo$algtype = "nonlinear"
  result$trfinfo  = trfinfo
  return(result)
}


# minimal addition --------------------------------------------------------
#' @keywords internal
#' @noRd
cge.minimaladd <- function(D0){
  ntgt = nrow(D0)
  if (as.integer(Matrix::rankMatrix(D0)) >= ntgt){
    output = list()
    output$D0 = D0
    output$alpha = 0.0
    return(output)
  } else {
    alpha = 0.1
    hello = TRUE
    while (hello){
      nice  = D0 + alpha*diag(ntgt)
      alpha = alpha + 0.1
      hello = (as.integer(Matrix::rankMatrix(nice)) < ntgt)
    }
    output = list()
    output$D0 = D0+alpha*diag(ntgt)
    output$alpha = alpha
    return(output)
  }
}

