#' Supervised Spectral Feature Selection
#'
#' SPEC algorithm selects features from the data via spectral graph approach.
#' Three types of ranking methods that appeared in the paper are available where
#' the graph laplacian is built via class label information.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of class labels.
#' @param ndim an integer-valued target dimension.
#' @param ranking types of feature scoring method. See the paper in the reference for more details.
#' @param preprocess an additional option for preprocessing the data. Default is "null". See also \code{\link{aux.preprocess}} for more details.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{sscore}{a length-\eqn{p} vector of spectral feature scores.}
#' \item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @examples
#' \donttest{
#' ## use iris data
#' ## it is known that feature 3 and 4 are more important.
#' data(iris)
#' set.seed(100)
#' subid = sample(1:150, 50)
#' iris.dat = as.matrix(iris[subid,1:4])
#' iris.lab = as.factor(iris[subid,5])
#'
#' ## try different ranking methods
#' out1 = do.specs(iris.dat, iris.lab, ranking="method1")
#' out2 = do.specs(iris.dat, iris.lab, ranking="method2")
#' out3 = do.specs(iris.dat, iris.lab, ranking="method3")
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, pch=19, col=iris.lab, main="SPECS::method1")
#' plot(out2$Y, pch=19, col=iris.lab, main="SPECS::method2")
#' plot(out3$Y, pch=19, col=iris.lab, main="SPECS::method3")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{zhao_spectral_2007}{Rdimtools}
#'
#' @seealso \code{\link{do.specu}}
#' @rdname feature_SPECS
#' @author Kisung You
#' @concept feature_methods
#' @export
do.specs <- function(X, label, ndim=2, ranking=c("method1","method2","method3"),
                     preprocess=c("null","center","scale","cscale","whiten","decorrelate")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label vector
  label  = check_label(label, n)
  #   3. ndim
  ndim = round(ndim)
  if (!check_ndim(ndim,p)){
    stop("* do.specs : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  #   4. preprocess
  if (missing(preprocess)){
    algpreprocess = "null"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   5. other parameters
  myscore = match.arg(ranking)

  #------------------------------------------------------------------------
  ## COMPUTATION : Preliminary
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #------------------------------------------------------------------------
  ## COMPUTATION : Graph Construction
  index = list()
  for (i in 1:length(unique(label))){
    index[[i]] = which(label==i)
  }
  indlength = rep(0,length(index))
  for (i in 1:length(index)){
    indlength[i] = length(index[[i]])
  }
  # affinity
  S = array(0,c(n,n))
  for (i in 1:(n-1)){
    labi = label[i]
    for (j in ((i+1):n)){
      labj = label[j]
      if (labi==labj){
        S[i,j] <- S[j,i] <- 1/indlength[labi]
      }
    }
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : normalized graph laplacian
  # normalized graph laplacian
  rowsumS  = base::rowSums(S)
  Dhalfvec = sqrt(rowsumS)
  Dhalfinv = diag(1/sqrt(rowsumS))
  L = diag(n) - Dhalfinv%*%S%*%Dhalfinv

  #------------------------------------------------------------------------
  ## COMPUTATION : case branching
  rankvec = switch(myscore,
                   "method1" = aux_spec_method1(pX, L, Dhalfvec),
                   "method2" = aux_spec_method2(pX, L, Dhalfvec),
                   "method3" = aux_spec_method3(pX, L, Dhalfvec))
  if (all(myscore=="method3")){
    idxvec = base::order(rankvec, decreasing=TRUE)[1:ndim]
  } else {
    idxvec = base::order(rankvec, decreasing=FALSE)[1:ndim]
  }
  projection = aux.featureindicator(p,ndim,idxvec)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$sscore  = rankvec
  result$featidx = idxvec
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}

# auxiliary functions -----------------------------------------------------
#' @keywords internal
aux_spec_method1 <- function(X, L, Dhalfvec){
  # preliminary computation
  p     = ncol(X)
  score = rep(0,p)

  # iteration
  for (i in 1:p){
    # compute normalized feature vector
    fi = Dhalfvec*as.vector(X[,i])
    fi = fi/sqrt(sum(fi^2))
    # compute the score
    score[i] = sum(as.vector(L%*%fi)*fi)
  }
  return(score)
}
#' @keywords internal
aux_spec_method2 <- function(X, L, Dhalfvec){
  # preliminary computation
  p     = ncol(X)
  score = rep(0,p)
  psi0  = as.vector(RSpectra::eigs(L, k=1, which="SR")$vectors) # smallest vector

  # iteration
  for (i in 1:p){
    # compute normalized feature vector
    fi = Dhalfvec*as.vector(X[,i])
    fi = fi/sqrt(sum(fi^2))
    # compute the terms
    term.top = sum(as.vector(L%*%fi)*fi)
    term.bot = 1 - sum(fi*psi0)
    # compute the score
    score[i] = term.top/term.bot
  }
  return(score)
}
#' @keywords internal
aux_spec_method3 <- function(X, L, Dhalfvec){
  # preliminary computation
  p     = ncol(X)
  score = rep(0,p)

  eigL    = base::eigen(L)
  eigLidx = (eigL$values > 100*.Machine$double.eps)

  evals   = eigL$values[eigLidx] # heuristic for (k-1) choice
  evecs   = eigL$vectors[,eigLidx]
  if (is.vector(evecs)){
    evecs = matrix(evecs, ncol=1)
  }
  nvals   = length(evals)

  # iteration
  for (i in 1:p){
    # compute normalized feature vector
    fi = Dhalfvec*as.vector(X[,i])
    fi = fi/sqrt(sum(fi^2))
    alpha = rep(0,nvals)
    for (j in 1:nvals){
      evecj    = as.vector(evecs[,j])
      term.top = sum(fi*evecj)
      term.bot = sqrt(sum(evecj^2))
      alpha[j] = term.top/term.bot
    }
    score[i] = sum((2-evals)*(alpha^2))
  }
  return(score)
}
