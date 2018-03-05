#' Double-Adjacency Graphs-based Discriminant Neighborhood Embedding
#'
#' Doublue Adjacency Graphs-based Discriminant Neighborhood Embedding (DAG-DNE) is a
#' variant of DNE. As its name suggests, it introduces two adjacency graphs for
#' homogeneous and heterogeneous samples accordaing to their labels.
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
#' ## generate data of 3 types with clear difference
#' dt1  = aux.gensamples(n=33)-100
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+100
#'
#' ## merge the data and create a label correspondingly
#' Y      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## try different numbers for neighborhood size
#' out1 = do.dagdne(Y, label, numk=5)
#' out2 = do.dagdne(Y, label, numk=10)
#' out3 = do.dagdne(Y, label, numk=25)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="nbd size=5")
#' plot(out2$Y[,1], out2$Y[,2], main="nbd size=10")
#' plot(out3$Y[,1], out3$Y[,2], main="nbd size=25")
#'
#' @references
#' \insertRef{ding_double_2015}{Rdimtools}
#'
#' @seealso \code{\link{do.dne}}
#' @author Kisung You
#' @rdname linear_DAGDNE
#' @export
do.dagdne <- function(X, label, ndim=2, numk=max(ceiling(nrow(X)/10),2),
                   preprocess=c("center","scale","cscale","decorrelate","whiten")){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. data matrix
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  #   2. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.dagdne : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   3. numk
  numk = as.integer(numk)
  if (!check_NumMM(numk,1,n/2,compact=FALSE)){stop("* do.dagdne : 'numk' should be an integer in [2,nrow(X)/2).")}
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
  #   2. Homogeneous and Heterogeneous Neighbors
  #   2-1. find neighborhood logical matrix
  logicalmat = aux.nbdlogical(pX, label, numk, numk)
  mat_hom    = logicalmat$hom
  mat_het    = logicalmat$het
  #   2-2. find adjacency matrix
  Fw = array(0,c(n,n))
  Fb = array(0,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      # 2-2-1. within-class adjacency matrix
      if ((mat_hom[i,j]==TRUE)||(mat_hom[j,i]==TRUE)){
        Fw[i,j] = 1
        Fw[j,i] = 1
      }
      # 2-2-2. intra-class adjacency matrix
      if ((mat_het[i,j]==TRUE)||(mat_het[j,i]==TRUE)){
        Fb[i,j] = 1
        Fb[j,i] = 1
      }
    }
  }
  #   2-3. build output matrix
  matS = (diag(rowSums(Fb))-Fb)-(diag(rowSums(Fw))-Fw)

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN COMPUTATION OF DAGDNE
  #   1. cost function
  costS = (t(pX)%*%matS%*%pX)
  #   2. RSpectra : largest
  projection = aux.adjprojection(RSpectra::eigs(costS, ndim)$vectors)


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
