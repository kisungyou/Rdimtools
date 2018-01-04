#' Semi-Supervised Discriminant Analysis
#'
#'
#'
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
#' #
#'
#' @author Kisung You
#' @export
do.sda <- function(X, label, ndim=2, type=c("proportion",0.1), alpha=1.0, beta=1.0){
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
  for (i in 1:length(ulabel)){
    if (sum(label==ulabel[i])==1){
      stop("* do.sda : no degerate class of size 1 is allowed.")
    }
  }
  if (all(!is.na(ulabel))){
    warning("* Semi-Supervised Learning : there is no missing labels. Consider using Supervised methods.")
  }
  labelorder = order(label)
  #   3. ndim
  if (!check_ndim(ndim,p)){
    stop("* do.sda : 'ndim' is a positive integer in [1,#(covariates)].")
  }
  ndim = as.integer(ndim)
  #   4. alpha : balancing
  alpha = as.double(alpha)
  if (!check_NumMM(alpha,0,Inf,compact=FALSE)){stop("* do.sda : 'alpha' needs to be a positive real number.")}
  #   5. beta : regularization
  beta = as.double(beta)
  if (!check_NumMM(beta,0,Inf,compact=TRUE)){stop("* do.sda : 'beta; needs to be a nonnegative real number.")}
  #   6. neighborhood type
  nbdtype = type

  #   (Implicit) preprocessing
  algpreprocess = "center"
  #   (Implicit) neighborhood symmetric
  nbdsymmetric = "union"

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. preprocessing with re-labeling of data
  tmplist = aux.preprocess(X,type=algpreprocess)
  trfinfo = tmplist$info
  pXoriginal   = tmplist$pX
  pX      = pXoriginal[labelorder,]
  trfinfo$algtype = "linear"
  label   = label[labelorder]
  #   2. neighborhood graph
  nbdstruct = aux.graphnbd(pX,method="euclidean",
                           type=nbdtype,symmetric=nbdsymmetric)
  nbdmask   = nbdstruct$mask
  #   3. S : binary adjacency
  S = nbdmask*1.0
  L = diag(rowSums(S))-S
  #   4. W : Weight Matrix
  idxmaxlabeled = sum(!is.na(label))
  Wl = sda_build_Wl(label[1:idxmaxlabeled])
  W  = array(0,c(n,n))
  W[1:idxmaxlabeled, 1:idxmaxlabeled] = Wl
  #   5. Itilde
  Itilde = array(0,c(n,n))
  Itilde[1:idxmaxlabeled, 1:idxmaxlabeled] = diag(idxmaxlabeled)

  #------------------------------------------------------------------------
  ## COMPUTATION : PRELIMINARY
  #   1. setup
  LHS = t(pX)%*%W%*%pX
  RHS = t(pX)%*%(    Itilde + (alpha*L) + (beta*diag(n)) )%*%pX
  #   2. top eigenvectors
  solgeigen = geigen::geigen(LHS, RHS, TRUE)$vectors[,p:(p-ndim+1)]
  #   3. adjust
  projection = aux.adjprojection(solgeigen)

  #------------------------------------------------------------------------
  ## RETURN
  result = list()
  result$Y = pXoriginal%*%projection
  result$trfinfo = trfinfo
  result$projection = projection
  return(result)





}




#  ------------------------------------------------------------------------
#' @keywords internal
#' @noRd
sda_build_Wl <- function(vec){
  uvec = unique(vec)
  c = length(uvec)
  l = length(vec)
  output = array(0,c(l,l))
  start = 1
  for (i in 1:c){
    li = sum(vec==uvec[i])
    output[start:(start+li-1),start:(start+li-1)] = 1/li
    start = start + li
  }
  return(output)
}







