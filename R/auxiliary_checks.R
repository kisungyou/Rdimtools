# LIST OF CHECKER FUNCTIONS
# 01. check_ndim   : universal
# 02. check_WMAT   : REE
# 03. check_label  : label information and return a well-processed one.
# 04. check_NumMM  : bounded number


# 01. check_ndim ----------------------------------------------------------
#' @noRd
#' @keywords internal
check_ndim <- function(ndim, p){
  if ((length(ndim)!=1)||(!is.numeric(ndim))||(ndim<1)||(ndim>p)||is.infinite(ndim)||is.na(ndim)){
    return(FALSE)
  } else {
    return(TRUE)
  }
}


# 02. check_WMAT ----------------------------------------------------------
#' @noRd
#' @keywords internal
check_WMAT <- function(W, n){
  # 1. size argument
  cond1 = ((is.matrix(W))&&(nrow(W)==n)&&(ncol(W)==n))
  # 2. no negative values
  cond2 = (all(W>=0))
  # 3. no Inf of NA
  cond3 = ((!any(is.na(W)))&&(!any(is.infinite(W))))

  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}



# 03. check_label ---------------------------------------------------------
# label should be corresponding to the number of observations.
#' @noRd
#' @keywords internal
check_label <- function(label, n){
  # 1. check if it is a proper vector
  if ((!is.vector(as.double(label)))||(length(label)!=n)){
    stop("* Supervised Learning : 'label' is required to be a vector of class labels.")
  }
  # 2. de-factoring the label
  label  = as.numeric(as.factor(label))
  ulabel = unique(label)
  K      = length(ulabel)
  if (K==1){
    warning("* Supervised Learning : 'label' should have at least 2 unique labelings.")
  }
  if (K==n){
    warning("* Supervised Learning : given 'label' has all unique elements.")
  }
  return(label)
}

# 04. check_NumMM ---------------------------------------------------------
# return TRUE if it is a valid one.
#' @noRd
#' @keywords internal
check_NumMM <- function(x, min, max, compact=TRUE){
  cond1 = (length(as.vector(x))==1)
  cond2 = (!is.na(x))
  if (compact){
    cond3 = ((min<=x)&&(x<=max))
  } else {
    cond3 = ((min<x)&&(x<max))
  }
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
