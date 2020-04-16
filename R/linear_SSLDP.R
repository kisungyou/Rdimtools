#' Semi-Supervised Locally Discriminant Projection
#'
#' Semi-Supervised Locally Discriminant Projection (SSLDP) is a semi-supervised
#' extension of LDP. It utilizes unlabeled data to overcome the small-sample-size problem
#' under the situation where labeled data have the small number. Using two information,
#' it both constructs the within- and between-class weight matrices incorporating the
#' neighborhood information of the data set.
#'
#' @examples
#' ## use iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.integer(iris$Species)
#'
#' ## copy a label and let 10% of elements be missing
#' nlabel = length(label)
#' nmissing = round(nlabel*0.10)
#' label_missing = label
#' label_missing[sample(1:nlabel, nmissing)]=NA
#'
#' ## compute with 3 different levels of 'beta' values
#' out1 = do.ssldp(X, label_missing, beta=0.1)
#' out2 = do.ssldp(X, label_missing, beta=0.5)
#' out3 = do.ssldp(X, label_missing, beta=0.9)
#'
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(out1$Y, col=label, main="SSLDP::beta=0.1")
#' plot(out2$Y, col=label, main="SSLDP::beta=0.5")
#' plot(out3$Y, col=label, main="SSLDP::beta=0.9")
#' par(opar)
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param label a length-\eqn{n} vector of data class labels.
#' @param ndim an integer-valued target dimension.
#' @param type a vector of neighborhood graph construction. Following types are supported;
#'  \code{c("knn",k)}, \code{c("enn",radius)}, and \code{c("proportion",ratio)}.
#'  Default is \code{c("proportion",0.1)}, connecting about 1/10 of nearest data points
#'  among all data points. See also \code{\link{aux.graphnbd}} for more details.
#' @param preprocess an additional option for preprocessing the data.
#' Default is "center". See also \code{\link{aux.preprocess}} for more details.
#' @param beta balancing parameter for intra- and inter-class information in \eqn{[0,1]}.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' \item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
#' }
#'
#' @references
#' \insertRef{zhang_semisupervised_2011}{Rdimtools}
#'
#' @author Kisung You
#' @rdname linear_SSLDP
#' @concept linear_methods 
#' @export
do.ssldp <- function(X, label, ndim=2, type=c("proportion",0.1),
                     preprocess=c("center","scale","cscale","whiten","decorrelate"), beta=0.5){
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
  if (any(is.infinite(ulabel))){
    stop("* Semi-Supervised Learning : no label of Inf is allowed.")
  }
  #   3. ndim
  ndim = as.integer(ndim)
  if (!check_ndim(ndim,p)){stop("* do.ssldp : 'ndim' is a positive integer in [1,#(covariates)).")}
  #   4. neighborhood information : asymmetric
  nbdtype = type
  nbdsymmetric = "union"
  #   5. preprocess
  if (missing(preprocess)){
    algpreprocess = "center"
  } else {
    algpreprocess = match.arg(preprocess)
  }
  #   6. beta
  beta = as.double(beta)
  if (!check_NumMM(beta,0,1,compact=TRUE)){stop("* do.ssldp : 'beta' is a balancing parameter in [0,1].")}

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing of data : note that output pX still has (n-by-p) format
  tmplist = aux.preprocess.hidden(X,type=algpreprocess,algtype="linear")
  trfinfo = tmplist$info
  pX      = tmplist$pX
  #   2. build neighborhood information
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask
  #   3. compute norm for squared difference norms and one's norms
  matDsq = (as.matrix(dist(pX))^2)
  vecnom = rep(0,n)
  for (i in 1:n){
    thevector = as.vector(pX[i,])
    vecnom[i] = sqrt(sum(thevector*thevector))
  }
  #   4. number of classes
  C = length(setdiff(ulabel,NA))

  #------------------------------------------------------------------------
  ## COMPUTATION : MAIN PART FOR SEMI-SUPERVISED LDP
  #   1. build weight matrices
  Hw = array(0,c(n,n))
  Hb = array(0,c(n,n))
  for (i in 1:(n-1)){
    class1 = label[i]
    for (j in (i+1):n){
      class2 = label[j]
      #   1-1. within-class weight
      if ((nbdmask[i,j]==TRUE)||(nbdmask[j,i]==TRUE)){
        if (((!is.na(class1))&&(!is.na(class2)))&&(class1==class2)){
          thevalue = exp(-((matDsq[i,j])/(vecnom[i]*vecnom[j])))
          Hw[i,j]  = thevalue
          Hw[j,i]  = thevalue
        } else if (is.na(class1)||is.na(class2)){
          thevalue = exp(-((matDsq[i,j])/(vecnom[i]*vecnom[j])))/C
          Hw[i,j]  = thevalue
          Hw[j,i]  = thevalue
        } else if (is.na(class1)&&is.na(class2)){
          thevalue = exp(-((matDsq[i,j])/(vecnom[i]*vecnom[j])))/(C^2)
          Hw[i,j]  = thevalue
          Hw[j,i]  = thevalue
        }
        #   1-2. between-class weight
        if (((!is.na(class1))&&(!is.na(class2)))&&(class1!=class2)){
          thevalue = exp(-((matDsq[i,j])/(vecnom[i]*vecnom[j])))
          Hb[i,j]  = thevalue
          Hb[j,i]  = thevalue
        }
      }
    }
  }
  #   2. build auxiliary matrices
  Dw = diag(rowSums(Hw))
  Lb = diag(rowSums(Hb))-Hb
  #   3. build cost functions
  LHS = t(pX)%*%(Lb-beta*Hw)%*%pX
  RHS = t(pX)%*%Dw%*%pX
  #   4. find projection : use top eigenvectors for geigen problem
  projection = aux.geigen(LHS, RHS, ndim, maximal=TRUE)


  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pX%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)
}
