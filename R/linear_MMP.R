#' Maximum Margin Projection
#'
#' Maximum Margin Projection (MMP) is a supervised linear method that maximizes the margin
#' between positive and negative examples at each local neighborhood based on
#' same- and different-class neighborhoods depending on class labels.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param numk the number of neighboring points.
#' @param alpha balancing parameter in \eqn{[0,1]}.
#' @param gamma weight for same-label data points with large magnitude.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=20)-100
#' dt2  = aux.gensamples(n=20)
#' dt3  = aux.gensamples(n=20)+100
#'
#' ## merge the data and create a label correspondingly
#' X      = rbind(dt1,dt2,dt3)
#' label  = rep(1:3, each=20)
#'
#' ## copy a label and let 20% of elements be missing
#' nlabel = length(label)
#' nmissing = round(nlabel*0.20)
#' label_missing = label
#' label_missing[sample(1:nlabel, nmissing)]=NA
#'
#' ## compare with PCA case for full-label case
#' ## for missing label case from MMP computation
#' out1 = do.pca(X, ndim=2)
#' out2 = do.mmp(X, label_missing, numk=10)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' plot(out1$Y, col=label, main="PCA projection")
#' plot(out2$Y, col=label, main="20% missing labels")
#' par(opar)
#'
#' @references
#' \insertRef{xiaofeihe_learning_2008}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_MMP
#' @concept linear_methods
#' @export
do.mmp <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"),
                   numk=max(ceiling(nrow(X)/10),2), alpha=0.5, gamma=50){
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
  if (all(!is.na(ulabel))){
    message("* Semi-Supervised Learning : there is no missing labels. Consider using Supervised methods.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){    stop("* do.mmp : 'ndim' is a positive integer in [1,#(covariates)].")  }
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }
  #   5. numk
  numk = as.integer(numk)
  if (!check_NumMM(numk,1,n/2,compact=FALSE)){stop("* do.mmp : 'numk' should be an integer in [2,nrow(X)/2).")}
  #   6. alpha
  alpha = as.double(alpha)
  if (!check_NumMM(alpha,0,1,compact=TRUE)){stop("* do.mmp : 'alpha' should be in [0,1].")}
  #   7. gamma
  gamma = as.double(gamma)
  if (!check_NumMM(gamma,0,1e+10,compact=FALSE)){stop("* do.mmp : 'gamma' should be a large positive real number.")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. perform PCA
  eigtest = eigen(cov(pX), only.values=TRUE)
  pcadim  = sum(eigtest$values > 0)
  if (pcadim <= ndim){
    warning("* do.mmp : target 'ndim' is larger than intrinsic data dimension achieved from PCA.")
    projection_first = diag(p)
    pcapX = pX
  } else{
    projection_first = aux.adjprojection(eigen(cov(pX))$vectors[,1:pcadim])
    pcapX = pX%*%projection_first
  }
  #   3. find neighborhood information with asymmetry
  nbdtype   = c("knn",numk)
  nbdstruct = aux.graphnbd(pcapX,method="euclidean",
                           type=nbdtype,symmetric="asymmetric")
  nbdmask   = nbdstruct$mask
  #   4. find Nb and Nw
  Nb = array(FALSE,c(n,n))
  Nw = array(FALSE,c(n,n))
  for (i in 1:n){
    # 4-1. current index set
    currentnbd = which(as.vector(nbdmask[i,]))
    currentlabels = label[currentnbd]
    # 4-2. Nb are neighbors with different labels
    idxtmp = which(currentlabels!=label[i])
    Nb[i,currentnbd[idxtmp]] = TRUE
    # 4-3. Nw are rest of it
    Nw[i,currentnbd[-idxtmp]] = TRUE
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR MMP
  #   1. build adjacency matrices
  #   1-1. Wb : neighbors with different labels
  Wb = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if ((Nb[i,j]==TRUE)||(Nb[j,i]==TRUE)){
        Wb[i,j] = 1.0
        Wb[j,i] = 1.0
      }
    }
  }
  #   1-2. Ww : rest of the diff-label points
  Ww = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (((!is.na(label[i]))&&(!is.na(label[j])))&&(label[i]==label[j])){
        Ww[i,j] = gamma
        Ww[j,i] = gamma
      } else if ((is.na(label[i])||is.na(label[j]))&&((Nw[i,j]==TRUE)||(Nw[j,i]==TRUE))){
        Ww[i,j] = 1.0
        Ww[j,i] = 1.0
      }
    }
  }
  #   2. build laplacian and others
  Db = diag(rowSums(Wb))
  Lb = Db - Wb
  Dw = diag(rowSums(Ww))
  #   3. cost function
  LHS = t(pcapX)%*%(alpha*Lb + (1-alpha)*Ww)%*%pcapX
  RHS = t(pcapX)%*%(Dw)%*%pcapX
  #   4. solve G-EVD with maximal solutions
  projection_second = aux.geigen(LHS, RHS, ndim)

  #------------------------------------------------------------------------
  ## RETURN
  #   1. merge two projection matrix
  projection = (projection_first%*%projection_second)

  #   2. return output
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
