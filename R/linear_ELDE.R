#' Exponential Local Discriminant Embedding
#'
#' Local Discriminant Embedding (LDE) suffers from a small-sample-size problem where
#' scatter matrix may suffer from rank deficiency. Exponential LDE (ELDE) provides
#' not only a remedy for the problem using matrix exponential, but also a flexible
#' framework to transform original data into a new space via distance diffusion mapping
#' similar to kernel-based nonlinear mapping.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param t kernel bandwidth in \eqn{(0,\infty)}.
#' @param preprocess  an additional option for preprocessing the data.
#' Default is "center" and other options of "decorrelate" and "whiten"
#' are supported. See also \code{\link{aux.preprocess}} for more details.
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
#' \dontrun{
#' ## generate data of 3 types with difference
#' dt1  = aux.gensamples(n=33)-100
#' dt2  = aux.gensamples(n=33)
#' dt3  = aux.gensamples(n=33)+100
#'
#' ## merge the data and create a label correspondingly
#' X      = rbind(dt1,dt2,dt3)
#' label  = c(rep(1,33), rep(2,33), rep(3,33))
#'
#' ## try different kernel bandwidth
#' out1 = do.elde(X, label, t=1)
#' out2 = do.elde(X, label, t=10)
#' out3 = do.elde(X, label, t=100)
#'
#' ## visualize
#' par(mfrow=c(1,3))
#' plot(out1$Y[,1], out1$Y[,2], main="ELDE::bandwidth=1")
#' plot(out2$Y[,1], out2$Y[,2], main="ELDE::bandwidth=10")
#' plot(out3$Y[,1], out3$Y[,2], main="ELDE::bandwidth=100")
#' }
#'
#' @references
#' \insertRef{dornaika_exponential_2013}{Rdimtools}
#'
#' @seealso \code{\link{do.lde}}
#' @rdname linear_ELDE
#' @author Kisung You
#' @export
do.elde <- function(X, label, ndim=2, t=1.0, preprocess=c("center","decorrelate","whiten"),
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
  if (!check_ndim(ndim,p)){stop("* do.elde : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. t
  t = as.double(t)
  if (!check_NumMM(t,.Machine$double.eps,Inf,compact=TRUE)){stop("* do.elde : 't' should be a positive real number.")}
  #   5. preprocess
  if (missing(preprocess)){    algpreprocess = "center"  }
  else {    algpreprocess = match.arg(preprocess)  }
  #   6. k1 and k2
  k1 = as.integer(k1)
  k2 = as.integer(k2)
  if (!check_NumMM(k1,1,n/2,compact=FALSE)){stop("* do.elde : 'k1' should be an integer in [2,nrow(X)/2).")}
  if (!check_NumMM(k2,1,n/2,compact=FALSE)){stop("* do.elde : 'k2' should be an integer in [2,nrow(X)/2).")}


  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pX      = tmplist$pX
  trfinfo$algtype = "linear"
  #   2. neighborhood information
  logicalmat = aux.nbdlogical(pX, label, k1, k2)
  Gw = logicalmat$hom
  Gb = logicalmat$het
  #   3. computing exp() distance
  expD = exp(-(as.matrix(dist(pX))^2)/t)

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR E-LDE
  #   1. Weight Matrices W
  Ww = expD*Gw
  Wb = Gb*1.0
  #Wb = expD*Gb
  #   2. compute auxiliary matrices
  Lw = diag(rowSums(Ww))-Ww
  Lb = diag(rowSums(Wb))-Wb

  Sw = t(pX)%*%Lw%*%pX
  Sb = t(pX)%*%Lb%*%pX

  #   3. scaling of the matrix
  Sw = Sw/Matrix::norm(Sw,"f")
  Sb = Sb/Matrix::norm(Sb,"f")

  #   3. matrix exponential
  expSw = as.matrix(Matrix::expm(Sw))
  expSb = as.matrix(Matrix::expm(Sb))


  #   4. compute projection matrix
  projection = aux.geigen(expSb, expSw, ndim, maximal=TRUE)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
