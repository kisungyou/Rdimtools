# LIST OF AUXILIARY FUNCTIONS
# 00. aux.typecheck        : check whether the data is poorly given
# 01. aux.preprocess -----------------------------------------------------------------------------
# 02. aux.gensamples -----------------------------------------------------------------------------
# 03. aux.graphnbd   -----------------------------------------------------------------------------
#     aux.graphnbdD        : distance matrix is already given
# 04. aux.shortestpath ---------------------------------------------------------------------------
# 05. aux.MaxMinLandmark   : choose a single landmark point
# 06. aux.kernelcov ------------------------------------------------------------------------------
# 07. aux.eigendec         : use Armadillo in a descending order
# 08. aux.pkgstat --------------------------------------------------------------------------------
# 09. aux.kernelcentering  : centering the kernel/gram matrix
# 10. aux.kernelprojection : given uncentered gram matrix, find the projected data
#                            note that it results (ndim-by-N) matrix, columns are projected vectors.
# 11. aux.adjprojection    : adjust projection matrix by simply normalizing each column
#     aux.adjqr            : qr decomposition is used for adjusting
# 12. aux.nbdlogical       : find homogeneous and heterogeneous neighborhood indexing
# 13. aux.geigen           : geigen in my taste
# 14. aux.featureindicator : generate (p-by-ndim) indicator matrix for projection
# 15. aux.traceratio.max   : compute trace ratio problem for maximal basis
# 16. aux.pinv             : use SVD and NumPy scheme
# 17. aux.bicgstab         : due to my stupidity, now Rlinsolve can't be used.
# 18. aux.oosprocess       : data processign for oos prediction
# 19. aux.findmaxidx       : find the row and column index of maximal elements
# 20. aux.randpartition    : given 1:n, divide it into K random partitions without replacement
# 21. aux.which.mink       : returns index of smallest
#     aux.which.maxk
# 22. aux.traceratio       : solve trace ratio problem with 2012 Ngo's algorithm
# 23. aux.2scatter         : compute LDA type within- and between- scatter matrices
# 24. aux.subsetid         : generate 'k' subset id

#  ------------------------------------------------------------------------
# 0. AUX.TYPECHECK
#  ------------------------------------------------------------------------
#' @noRd
#' @keywords internal
aux.typecheck <- function(data, verbose=FALSE){
  # data frame into matrix
  matinput = as.matrix(data)
  if ((dim(matinput)[1]==1)||(dim(matinput)[2]==1)){
    if (verbose==TRUE){
      message("WARNING : input data should be matrix, not a vector.")
    }
    matinput = matrix(matinput)
  }

  # check if there exists any NA, Inf, -Inf
  if (!is.numeric(matinput)){
    warning("ERROR : input data should be numeric arrays.")
    return(FALSE)
  }
  if (any(is.infinite(matinput))||any(is.na(matinput))){
    warning("ERROR : input data should contain no Inf or NA values.")
    return(FALSE)
  }
  return(TRUE)
}



#  ------------------------------------------------------------------------
# 5. AUX.MAXMINLANDMARK
#  ------------------------------------------------------------------------
#' @noRd
#' @keywords internal
aux.MaxMinLandmark <- function(X,npoints,pdflag=FALSE){
  # 5-1. setting
  npoints = as.integer(npoints)
  if (nrow(X)<=npoints){
    stop("ERROR : npoints should be smaller than the total number of original data points.")
  }

  # 5-2. initialize with two starting points
  landmarkidx = array(0,c(1,npoints))
  nX = nrow(X)
  seqnp = seq_len(nX)
  idx12 = sample(1:nX,2)

  landmarkidx[1:2] = idx12     # recorded random points
  seqnp = setdiff(seqnp,idx12) # sequence left over

  # 5-3. main computation
  if (pdflag==TRUE){
    pD = X
  } else {
    pD = as.matrix(dist(X))
  }
  if (npoints>2){
    for (it in 3:npoints){
      # 5-3-1. marginalize the data
      plandmark  = landmarkidx[1:(it-1)]
      currentidx = aux_landmarkMaxMin(pD, plandmark, seqnp)

      # 5-3-2. update
      landmarkidx[it] = currentidx
      seqnp = setdiff(seqnp,currentidx)
    }
  }

  # 5-4. return results
  return(landmarkidx)
}


# 7. eigendecomposition : Armadillo + Descending order --------------------
#' @noRd
#' @keywords internal
aux.eigendec <- function(X){
  if (nrow(X)!=ncol(X)){
    stop("ERROR : a given matrix is not a square form.")
  }
  output = aux_eigendecomposition(X);

  result = list()
  result$eigval = rev(output$eigval)
  result$eigvec = output$eigvec[,rev(seq_len(length(output$eigval)))]
  return(result)
}

# 09. aux.kernelcentering -------------------------------------------------
#     centering the kernel/gram matrix
#' @keywords internal
#' @noRd
aux.kernelcentering <- function(K){
  N = nrow(K)
  if (ncol(K)!=N){
    stop("* aux.kernelcentering : an input K should be a square matrix.")
  }
  onesN  = array(1,c(N,N))/N
  Ktilde = K-(onesN%*%K)-(K%*%onesN)+(onesN%*%K%*%onesN)
  return(Ktilde)
}


# 10. aux.kernelprojection ------------------------------------------------
# given uncentered gram matrix, find the projected data
# note that it returns (ndim-by-N)
#' @keywords internal
#' @noRd
aux.kernelprojection <- function(KK, ndim){
  KKcentered = aux.kernelcentering(KK)
  KKceigen   = eigen(KKcentered)
  Y          = (t(KKceigen$vectors[,1:ndim]) %*% KK)
  return(Y)
}


# 11. aux.adjprojection : adjust projection matrix ------------------------
#' @keywords internal
#' @noRd
aux.adjqr <- function(P){
  p = ncol(P)
  Pid = (t(P)%*%P)
  dnormval = base::norm(diag(p)-Pid, type = "F")
  if (dnormval > 1e-10){
    return(qr.Q(qr(P)))
  } else {
    return(P)
  }
}
#' @keywords internal
#' @noRd
aux.adjprojection <- function(P){
  PP = aux.adjqr(P)
  n = nrow(P)
  p = ncol(P)
  output = array(0,c(n,p))
  for (i in 1:p){
    cvec = as.vector(P[,i])
    output[,i] = cvec/sqrt(sum(cvec*cvec))
  }
  return(output)
}
# 12. aux.nbdlogical : find homogeneous and heterogeneous neighbor --------
#' @keywords internal
#' @noRd
aux.nbdlogical <- function(X, label, khomo, khet){
  D = as.matrix(dist(X))
  n = nrow(D)
  # 1. homogeneous logical matrix
  logical_hom = array(FALSE,c(n,n))
  for (i in 1:n){
    # 1-1. index of same class
    idxhom = setdiff(which(label==label[i]), i)
    # 1-2. which is the smallest numk ?
    partD = as.vector(D[i,idxhom])
    # 1-3. partially smallest ones
    partidx = which(      partD <= max(sort(partD)[1:max(min(khomo, length(idxhom)),1)])    )
    # 1-4. adjust idxhom
    idxhomadj = idxhom[partidx]
    logical_hom[i,idxhomadj] = TRUE
  }
  # 2. heterogeneous logical matrix
  logical_het = array(FALSE,c(n,n))
  for (i in 1:n){
    idxhet = which(label!=label[i])
    partD = as.vector(D[i, idxhet])
    partidx = which(partD <= max(sort(partD)[1:max(min(khet, length(idxhet)),1)]))
    idxhetadj = idxhet[partidx]
    logical_het[i,idxhetadj] = TRUE
  }
  output = list()
  output$hom = logical_hom
  output$het = logical_het
  return(output)
}

# 13. aux.geigen : now it uses RcppArmadillo ------------------------------
#' @keywords internal
#' @noRd
aux.geigen <- function(top, bottom, ndim, maximal=TRUE){
  # new 1. use 'maotai's mimic function : increasing order as well
  gfun  = getFromNamespace("hidden_geigen","maotai")
  geigs = gfun(top, bottom, normalize=TRUE)
  nobj  = length(geigs$values)
  ndim  = round(ndim)

  # new 2. separate eigenvectors accordingly to the orders
  if (maximal){
    vectors = geigs$vectors[,nobj:(nobj-ndim+1)]
  } else {
    vectors = geigs$vectors[,1:ndim]
  }

  # new 3. return
  if (ndim==1){
    return(matrix(vectors, ncol=1))
  } else {
    return(vectors)
  }

  # # 1. first, run CPP with RcppArmadillo -> change to 'geigen' package
  # # geigs = aux_geigen(top, bottom) # Armadillo goes Decreasing order
  # geigs = geigen::geigen(top, bottom) # geigen goes Increasing order
  # maxp  = length(geigs$values)
  #
  # # 2. separate values and vectors; change correspondingly for Decreasing order
  # values  = geigs$values[maxp:1]
  # vectors = geigs$vectors[,maxp:1]
  #
  # if (maximal==TRUE){
  #   partvals = values[1:ndim]
  # } else {
  #   partvals = values[maxp:(maxp-ndim+1)]
  # }
  # if (all(base::abs(base::Im(partvals))<(100*(.Machine$double.eps)))){
  #   values  = base::Re(values)
  #   vectors = aux.adjprojection(base::Re(vectors))
  # } else {
  #   stop("* aux.geigen : generalized eigenvalue problem returned imaginary eigenvalues.")
  # }
  #
  #
  # # 3. branching case
  # if (ndim > 1){
  #   if (maximal==TRUE){
  #     projection = vectors[,1:ndim]
  #   } else {
  #     projection = vectors[,maxp:(maxp-ndim+1)]
  #   }
  # } else {
  #   if (maximal==TRUE){
  #     vecsol = vectors[,1]
  #     projection = matrix(vecsol/sqrt(sum(vecsol*vecsol)))
  #   } else {
  #     vecsol = vectors[,maxp]
  #     projection = matrix(vecsol/sqrt(sum(vecsol*vecsol)))
  #   }
  # }
  # return(projection)
}


# 14. aux.featureindicator : generate indicator matrix --------------------
#     generate (p-by-ndim) indicator matrix for projection
#' @keywords internal
#' @noRd
aux.featureindicator <- function(p,ndim,idxvec){
  if (length(idxvec)!=ndim){
    stop("* aux.featureindicator : selection had some problem.")
  }
  output = array(0,c(p,ndim))
  for (i in 1:ndim){
    selectedcolumn = as.integer(idxvec[i])
    output[selectedcolumn,i] = 1
  }
  return(output)
}

# 15. aux.traceratio.max : compute trace ratio problem for maximal --------
#' @keywords internal
#' @noRd
aux.traceratio.max <- function(A,B,ndim,tol=1e-10){
  # -------------------------------------------------------------
  # PRELIMINARY SETTING
  if (nrow(A)!=ncol(A)){
    stop("* aux.traceratio : A is not square.")
  } else {
    d = nrow(A)
  }
  if (nrow(B)!=ncol(B)){
    stop("* aux.traceratio : B is not square")
  }
  if (nrow(B)!=d){
    stop("* aux.traceratio : B is not matching size of A.")
  }
  m = round(ndim)
  r = round(aux_rank(B)) # as.integer(Matrix::rankMatrix (B))

  # -------------------------------------------------------------
  # MAIN BRANCHING
  if (m <= (d-r)){
    Z = RSpectra::eigs(B,(d-r),which="SR")$vectors
    cost = t(Z)%*%A%*%Z
    Zright = RSpectra::eigs(cost, m)$vectors
    output = Z%*%Zright;
  } else {
    lambda1 = sum(diag(A))/sum(diag(B))
    lambda2 = sum(base::eigen(A, only.values=TRUE)$values[1:m])/sum(base::eigen(B, only.values = TRUE)$values[1:m])
    lambda  = ((lambda1+lambda2)/2)
    while ((lambda2 -lambda1)>tol){
      gamma = sum(base::eigen((A-lambda*B), only.values = TRUE)$values[1:m])
      if (gamma > 0){
        lambda1 = lambda
      } else {
        lambda2 = lambda
      }
      lambda = ((lambda1+lambda2)/2)
    }
    output = RSpectra::eigs((A-lambda*B),ndim)$vectors
  }
  return(output)
}

# 16. aux.pinv : use SVD and NumPy Scheme ---------------------------------
# https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Singular_value_decomposition_(SVD)
#' @keywords internal
#' @noRd
aux.pinv <- function(A){
  svdA      = base::svd(A)
  tolerance = (.Machine$double.eps)*max(c(nrow(A),ncol(A)))*as.double(max(svdA$d))

  idxcut    = which(svdA$d <= tolerance)
  invDvec   = (1/svdA$d)
  invDvec[idxcut] = 0

  output = (svdA$v%*%diag(invDvec)%*%t(svdA$u))
  return(output)
}




# 17. aux.bicgstab --------------------------------------------------------
#' @keywords internal
#' @noRd
aux.bicgstab <- function(A,B,xinit=NA,reltol=1e-5,maxiter=1000,
                            preconditioner=diag(ncol(A)),verbose=TRUE){
  ###########################################################################
  # Step 0. Initialization
  if (verbose){
    message("* aux.bicgstab : Initialiszed.")
  }
  if (any(is.na(A))||any(is.infinite(A))||any(is.na(B))||any(is.infinite(B))){
    stop("* aux.bicgstab : no NA or Inf values allowed.")
  }

  A = matrix(A,nrow=nrow(A))
  if (is.vector(B)){
    B = matrix(B,ncol=1)
  } else {
    B = matrix(B,nrow=nrow(B))
  }
  preconditioner = matrix(preconditioner,nrow=nrow(preconditioner))
  sparseflag = FALSE

  # sparseformats = c("dgCMatrix","dtCMatrix","dsCMatrix")
  # if ((class(A)%in%sparseformats)||(class(B)%in%sparseformats)||(class(preconditioner)%in%sparseformats)){
  #   A = Matrix(A,sparse=TRUE)
  #   B = Matrix(B,sparse=TRUE)
  #   preconditioner = Matrix(preconditioner,sparse=TRUE)
  #   sparseflag = TRUE
  # } else {
  #   A = matrix(A,nrow=nrow(A))
  #   if (is.vector(B)){
  #     B = matrix(B)
  #   } else {
  #     B = matrix(B,nrow=nrow(B))
  #   }
  #   preconditioner = matrix(preconditioner,nrow=nrow(preconditioner))
  #   sparseflag = FALSE
  # }
  # xinit
  if (is.na(xinit)){
    xinit = matrix(rnorm(ncol(A)))
  } else {
    if (length(xinit)!=ncol(A)){
      stop("* aux.bicgstab : 'xinit' has invalid size.")
    }
    xinit = matrix(xinit)
  }
  ###########################################################################
  # Step 1. Preprocessing
  # 1-1. Neither NA nor Inf allowed.
  if (any(is.infinite(A))||any(is.na(A))||any(is.infinite(B))||any(is.na(B))){
    stop("* aux.bicgstab : no NA, Inf, -Inf values are allowed.")
  }
  # 1-2. Size Argument
  m = nrow(A)
  if (is.vector(B)){
    mB = length(B)
    if (m!=mB){
      stop("* aux.bicgstab : a vector B should have a length of nrow(A).")
    }
  } else {
    mB = nrow(B)
    if (m!=mB){
      stop("* aux.bicgstab : an input matrix B should have the same number of rows from A.")
    }
  }
  if (is.vector(B)){
    B = as.matrix(B)
  }
  # 1-3. Adjusting Case
  if (m > ncol(A)){        ## Case 1. Overdetermined
    B = t(A)%*%B
    A = t(A)%*%A
  } else if (m < ncol(A)){ ## Case 2. Underdetermined
    stop("* aux.bicgstab : underdetermined case is not supported.")
  }
  # 1-4. Preconditioner : only valid for square case
  if (!all.equal(dim(A),dim(preconditioner))){
    stop("* aux.bicgstab : Preconditioner is a size-matching.")
  }
  if (verbose){message("* aux.bicgstab : preprocessing finished ...")}
  ###########################################################################
  # Step 2. Main Computation
  ncolB = ncol(B)
  if (ncolB==1){
    if (!sparseflag){
      vecB = as.vector(B)
      res = linsolve.bicgstab.single(A,vecB,xinit,reltol,maxiter,preconditioner)
    } else {
      vecB = B
      res = linsolve.bicgstab.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner)
    }
  } else {
    x      = array(0,c(ncol(A),ncolB))
    iter   = array(0,c(1,ncolB))
    errors = list()
    for (i in 1:ncolB){
      vecB = as.vector(B[,i])
      tmpres = linsolve.bicgstab.single(A,vecB,xinit,reltol,maxiter,preconditioner)
      # if (!sparseflag){
      #
      # } else {
      #   vecB = Matrix (B[,i],sparse=TRUE)
      #   tmpres = linsolve.bicgstab.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner)
      # }
      x[,i]        = tmpres$x
      iter[i]      = tmpres$iter
      errors[[i]]  = tmpres$errors
      if (verbose){
        message(paste("* aux.bicgstab : B's column.",i,"being processed.."))
      }
    }
    res = list("x"=x,"iter"=iter,"errors"=errors)
  }

  ###########################################################################
  # Step 3. Finalize
  if ("flag" %in% names(res)){
    flagval = res$flag
    if (flagval==0){
      if (verbose){
        message("* aux.bicgstab : convergence well achieved.")
      }
    } else if (flagval==1){
      if (verbose){
        message("* aux.bicgstab : convergence not achieved within maxiter.")
      }
    } else {
      if (verbose){
        message("* aux.bicgstab : breakdown.")
      }
    }
    res$flag = NULL
  }
  if (verbose){
    message("* aux.bicgstab : computations finished.")
  }
  return(res)
}

# 18. aux.oospreprocess ---------------------------------------------------
#     data processing for out-of-sample prediction
#' @keywords internal
#' @noRd
aux.oospreprocess <- function(data, trfinfo){
  ## 0. parameter
  n = nrow(data)
  p = ncol(data)
  output = array(0,c(n,p))
  ## 1. extract mean
  meanvec = as.vector(trfinfo$mean)
  for (i in 1:n){
    output[i,] = (as.vector(data[i,])-meanvec)
  }
  ## 2. multiplier
  multiplier = trfinfo$multiplier
  if (is.matrix(multiplier)){
    output = output%*%multiplier
  } else {
    output = output*multiplier
  }
  return(output)
}

# 19. aux_findmaxidx       : find the row and column index of maxi --------
#' @keywords internal
#' @noRd
aux.findmaxidx <- function(A){
  # 19. aux_findmaxidx       : find the row and column index of maximal elements
  output = which(A==max(A), arr.ind=TRUE)
  if (nrow(output)>1){
    return(output[1,])
  } else {
    return(output)
  }
}


# 20. aux.randpartition ---------------------------------------------------
#     given 1:n, divide it into K random partitions without replacement
#' @keywords internal
#' @noRd
aux.randpartition <- function(n, K){
  output = list()
  if (K==1){
    output = list()
    output[[1]] = 1:n
  } else {
    listall = 1:n
    singleK = round(n/K)
    for (i in 1:(K-1)){
      output[[i]] = sample(listall, singleK, replace = FALSE)
      listall = setdiff(listall, output[[i]])
    }
    output[[K]] = listall
  }
  return(output)
}


# 21. aux.which.mink ------------------------------------------------------
#     aux.which.maxk
#' @keywords internal
#' @noRd
aux.which.mink <- function(x, k=1){
  return(order(x)[1:k])
}
#' @keywords internal
#' @noRd
aux.which.maxk <- function(x, k=1){
  return(order(x,decreasing = TRUE)[1:k])
}


# 22. aux.traceratio  : solve trace ratio problem with 2012 Ngo's  --------
#' @keywords internal
#' @noRd
aux.traceratio <- function(A, B, dim, eps, maxiter){
  ## in the the language
  n = nrow(A)
  p = dim

  ## prepare the initializer
  Vold = qr.Q(qr(matrix(rnorm(n*p),ncol=p)))
  rhoold = 0
  for (i in 1:maxiter){
    Vnew   = RSpectra::eigs(A-rhoold*B,p,which="LR")$vectors
    rhonew = sum(diag(t(Vnew)%*%A%*%Vnew))/sum(diag(t(Vnew)%*%B%*%Vnew))

    rhoinc = abs(rhonew-rhoold)
    Vold   = Vnew
    rhoold = rhonew

    if (rhoinc < eps){
      break
    }
  }

  ## let's try to return !
  return(Vold)
}


# 23. aux.2scatter --------------------------------------------------------
#     data should be provided as a matrix (columns are variables)
#' @keywords internal
#' @noRd
aux.2scatter <- function(pX, label){
  # 0. extra information
  if (is.vector(pX)){
    n = length(pX)
    p = 1
  } else {
    n = nrow(pX)
    p = ncol(pX)
  }


  # 1. extract label information
  label   = round(label)
  ulabel  = unique(label)
  datlist = list()
  for (i in 1:length(ulabel)){
    if (is.vector(pX)){
      datlist[[i]] = pX[(label==ulabel[i])]
    } else {
      datlist[[i]] = pX[(label==ulabel[i]),]
    }
  }
  # 2. compute two types of scatter matrices
  #   2-1. error matrix/E/error variance
  scattermat <- function(x){
    return(cov(x)*(nrow(x)-1))
  }
  matE = array(0,c(p,p))
  for (i in 1:length(ulabel)){
    if (is.vector(datlist[[i]])){
      tgt = as.matrix(datlist[[i]])
    } else {
      tgt = datlist[[i]]
    }
    matE = matE + cov(tgt)*(nrow(tgt)-1)
  }
  matH = array(0,c(p,p))
  if (is.vector(datlist[[1]])){
    meanlist = lapply(datlist, base::mean)
  } else {
    meanlist = lapply(datlist, colMeans)
  }
  if (is.vector(pX)){
    meantott = base::mean(pX)
  } else {
    meantott = colMeans(pX)
  }
  for (i in 1:length(ulabel)){
    meandiff = as.vector(meanlist[[i]]-meantott)
    if (is.vector(datlist[[i]])){
      matH = matH + length(datlist[[i]])*outer(meandiff,meandiff)
    } else {
      matH = matH + nrow(datlist[[i]])*outer(meandiff,meandiff)
    }
  }

  output = list()
  if (is.vector(pX)){
    output$within = as.double(matE)
    output$between = as.double(matH)
  } else {
    output$within = matE
    output$between = matH
  }

  return(output)
}


# 24. aux.subsetid --------------------------------------------------------
#' @keywords internal
aux.subsetid <- function(n, k){
  x = sample(1:n)
  return(split(x, sort(x%%k)))
}
