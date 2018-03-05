#' Maximal Local Interclass Embedding
#'
#' Maximal Local Interclass Embedding (MLIE) is a linear supervised method that
#' the local interclass graph and the intrinsic graph are constructed to find a set of
#' projections that maximize the local interclass scatter and the local
#' intraclass compactness at the same time. It can be deemed an extended version of MFA.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param k1 the number of same-class neighboring points (homogeneous neighbors).
#' @param k2 the number of different-class neighboring points (heterogeneous neighbors).
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
#' X      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## try different numbers for neighborhood size
#' out1 = do.mlie(X, label, k1=5, k2=5)
#' out2 = do.mlie(X, label, k1=10,k2=10)
#' out3 = do.mlie(X, label, k1=25,k2=25)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="MLIE::nbd size=5")
#' plot(out2$Y[,1], out2$Y[,2], main="MLIE::nbd size=10")
#' plot(out3$Y[,1], out3$Y[,2], main="MLIE::nbd size=25")
#'
#' @references
#' \insertRef{lai_maximal_2011}{Rdimtools}
#'
#' @seealso \code{\link{do.mfa}}
#' @rdname linear_MLIE
#' @export
do.mlie <- function(X, label, ndim=2, preprocess=c("center","scale","cscale","decorrelate","whiten"),
                    k1=max(ceiling(nrow(X)/10),2), k2=max(ceiling(nrow(X)/10),2)){
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
  if (any(is.na(label))||(any(is.infinite(label)))){stop("* Supervised Learning : any element of 'label' as NA or Inf will simply be considered as a class, not missing entries.")  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.mlie : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }
  #   5. k1 and k2
  k1 = as.integer(k1)
  k2 = as.integer(k2)
  if (!check_NumMM(k1,1,n/2,compact=FALSE)){stop("* do.mlie : 'k1' should be an integer in [2,nrow(X)/2).")}
  if (!check_NumMM(k2,1,n/2,compact=FALSE)){stop("* do.mlie : 'k2' should be an integer in [2,nrow(X)/2).")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX

  #   2. PCA preprocessing
  eigtest = eigen(cov(pX), only.values=TRUE)
  pcadim  = sum(eigtest$values > 0)
  if (pcadim <= ndim){
    warning("* do.mlie : target 'ndim' is larger than intrinsic data dimension achieved from PCA.")
    projection_first = diag(p)
    pcapX = pX
  } else{
    projection_first = aux.adjprojection(RSpectra::eigs(cov(pX), pcadim)$vectors)
    pcapX = pX%*%projection_first
  }

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR MFA
  #   1. compute homogeneous (intraclass) and heterogeneous (interclass) neighborhood
  #   1-1. compute neighborhood
  logicalmat = aux.nbdlogical(pcapX, label, k1, k2)
  mat_hom    = logicalmat$hom
  mat_het    = logicalmat$het
  #   1-2. OR class under symmetrization
  W  = array(as.logical(mat_hom + t(mat_hom)),c(n,n))*1.0; diag(W)=0
  Ws = array(as.logical(mat_het + t(mat_het)),c(n,n))*1.0; diag(Ws)=0
  #   1-3. D and Ds as well
  Ds = diag(rowSums(Ws))
  D  = diag(rowSums(W))

  #   2. formulation as generalized eigenvalue problem
  LHS = t(pcapX)%*%(D-W)%*%pcapX
  RHS = t(pcapX)%*%(Ds-Ws)%*%pcapX

  #   3. use as MAXIMIZATION Problem
  projection_second = aux.geigen(LHS, RHS, ndim, maximal=TRUE)


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
