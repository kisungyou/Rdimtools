#' Supervised Laplacian Eigenmaps
#'
#' Supervised Laplacian Eigenmaps (SPLAPEIG) is a supervised variant of Laplacian Eigenmaps.
#' Instead of setting up explicit neighborhood, it utilizes an adaptive threshold strategy
#' to define neighbors for both within- and between-class neighborhood. It then builds affinity
#' matrices for each information and solves generalized eigenvalue problem. This algorithm
#' may be quite sensitive in the choice of \code{beta} value.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "null" and three options of "center", "decorrelate", or "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
#' @param beta bandwidth parameter for heat kernel in \eqn{[0,\infty)}.
#' @param gamma a balancing parameter in \eqn{[0,1]} between within- and between-class information.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=33)-100
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+100
#'
#' ## merge the data and create a label correspondingly
#' X      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## try different balancing parameters with beta=100
#' out1 = do.splapeig(X, label, beta=100, gamma=0.1)
#' out2 = do.splapeig(X, label, beta=100, gamma=0.5)
#' out3 = do.splapeig(X, label, beta=100, gamma=0.9)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="gamma=0.1")
#' plot(out2$Y[,1], out2$Y[,2], main="gamma=0.5")
#' plot(out3$Y[,1], out3$Y[,2], main="gamma=0.9")
#' }
#'
#'
#' @references
#' \insertRef{raducanu_supervised_2012}{Rdimtools}
#'
#' @seealso \code{\link{do.lapeig}}
#' @author Kisung You
#' @rdname nonlinear_SPLAPEIG
#' @export
do.splapeig <- function(X, label, ndim=2,
                        preprocess=c("null","center","whiten","decorrelate"),
                        beta=1.0, gamma=0.5){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. label : check and return a de-factored vector
  #   For this example, there should be no degenerate class of size 1.
  label  = check_label(label, n)
  if (any(is.na(label))||(any(is.infinite(label)))){stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.splapeig : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "null"  }
  else {    algpreprocess = match.arg(preprocess)  }
  #   5. beta : kernel parameter
  beta = as.double(beta)
  if (!check_NumMM(beta,0,Inf,compact=TRUE)){stop("* do.splapeig : 'beta' should be a nonnegative real number.")}
  #   6. gamma : balancing parameter
  gamma = as.double(gamma)
  if (!check_NumMM(gamma,0,1,compact=TRUE)){stop("* do.splapeig : 'gamma' is a balancing parameter in [0,1].")}


  #------------------------------------------------------------------------
  ## PREPROCESSING
  if (algpreprocess=="null"){
    trfinfo = list()
    trfinfo$type = "null"
    pX = as.matrix(X,nrow=nrow(X));
  } else {
    tmplist = aux.preprocess(X,type=algpreprocess)
    trfinfo = tmplist$info
    pX      = tmplist$pX
  }
  trfinfo$algtype = "nonlinear"

  #------------------------------------------------------------------------
  ## MAIN COMUTATION FOR SUPERVISED LAPLACIAN EIGENMAPS
  DmatSq  = (as.matrix(dist(pX))^2)
  expbeta = exp(-DmatSq/beta)
  #   1. compute adaptive threshold
  AS = rep(0,n)
  for (i in 1:n){
    #   1-1. select the vector and compute exp(-SQ/beta)
    rowvec = exp(-as.vector(DmatSq[i,])/beta)
    #   1-2. summation and division
    AS[i] = (sum(rowvec))/n
  }
  #   2. two neighborhoods
  Nw = array(FALSE, c(n,n))
  Nb = array(FALSE, c(n,n))
  for (i in 1:n){
    #   2-1. two labels
    idx_same = which(label==label[i])
    idx_diff = which(label!=label[i])
    #   2-2. larger than AS(y_i)
    cexpsq    = as.vector(expbeta[i,])
    idx_large = which(cexpsq > AS[i])
    #   2-3.
    idxW = setdiff(intersect(idx_same, idx_large),i)
    idxB = setdiff(intersect(idx_diff, idx_large),i)
    #   2-4. fill in the logicals
    Nw[i,idxW] = TRUE
    Nb[i,idxB] = TRUE
  }
  #   3. build affinity matrices
  Ww = array(0,c(n,n))
  Wb = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      #   3-1. Ww first; within-class
      if ((Nw[i,j]==TRUE)||(Nw[j,i]==TRUE)){
        Wwvalue = as.double(expbeta[i,j])
        Ww[i,j] = Wwvalue
        Ww[j,i] = Wwvalue
      }
      #   3-2. Wb next; between-class
      if ((Nb[i,j]==TRUE)||(Nb[j,i]==TRUE)){
        Wb[i,j] = 1.0
        Wb[j,i] = 1.0
      }
    }
  }
  #   4. build cost function
  #   4-1. materials
  Dw = diag(rowSums(Ww))
  Lb = diag(rowSums(Wb))-Wb
  #   4-2. B
  B = gamma*Lb + (1-gamma)*Ww

  #   5. compute TOP eigenvectors
  Youtput = aux.geigen(B, Dw, ndim, maximal=TRUE)



  #------------------------------------------------------------------------
  ## RETURN OUTPUT
  result = list()
  result$Y = Youtput
  result$trfinfo  = trfinfo
  return(result)
}


