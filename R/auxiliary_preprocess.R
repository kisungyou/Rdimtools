#' Preprocessing the data
#'
#' \code{aux.preprocess} can perform one of following operations; \code{"center"}, \code{"scale"},
#' \code{"cscale"}, \code{"decorrelate"} and \code{"whiten"}. See below for more details.
#'
#' @section
#' Operations: we have following operations,
#'  \describe{
#'  \item{\code{"center"}}{subtracts mean of each column so that every variable has mean \eqn{0}.}
#'  \item{\code{"scale"}}{turns each column corresponding to variable have variance \eqn{1}.}
#'  \item{\code{"cscale"}}{combines \code{"center"} and \code{"scale"}.}
#'  \item{\code{"decorrelate"}}{\code{"center"} and sets its covariance term having diagonal entries only.}
#'  \item{\code{"whiten"}}{\code{"decorrelate"} and sets all diagonal elements be \eqn{1}.}
#'  }
#'
#' @param data an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param type one of \code{"center"}, \code{"scale"}, \code{"cscale"}, \code{"decorrelate"} or \code{"whiten"}.
#' @return named list containing:
#' \describe{
#' \item{pX}{an \eqn{(n\times p)} matrix after preprocessing in accordance with \code{type} parameter}
#' \item{info}{a list containing \itemize{
#' \item \code{type:} name of preprocessing procedure.
#' \item \code{mean:} a mean vector of length \eqn{p}.
#' \item \code{multiplier:} a \eqn{(p\times p)} matrix or 1 for "center".}}
#' }
#'
#' @examples
#' \donttest{
#' ## Generate data
#' set.seed(100)
#' X = aux.gensamples(n=200)
#'
#' ## 5 types of preprocessing
#' X_center = aux.preprocess(X)
#' X_scale  = aux.preprocess(X,type="scale")
#' X_cscale = aux.preprocess(X,type="cscale")
#' X_decorr = aux.preprocess(X,type="decorrelate")
#' X_whiten = aux.preprocess(X,type="whiten")
#' }
#'
#' @rdname aux_preprocess
#' @author Kisung You
#' @export
aux.preprocess <- function(data,type=c("center","scale","cscale","decorrelate","whiten")){
  # data : (n-by-d)
  #   ARMA with CPP
  #   input & output are (d-by-n)
  #   Don't forget to transpose at the end.
  if (is.data.frame(data)){
    matinput = t(as.matrix(data))
  } else if (is.matrix(data)){
    matinput = t(data)
  } else {
    warning("WARNING : input should be either dataframe or matrix.")
  }
  type = match.arg(type)

  ## two added methods : 'scale' and 'cscale'
  if (type=="scale"){ # add 1 : "scale"
    tmpdata = t(matinput)
    p       = ncol(tmpdata)

    multiplier = rep(0,p)
    for (i in 1:p){
      multiplier[i] = 1/stats::sd(as.vector(tmpdata[,i]))
    }
    info = list()
    info$type = "scale"
    info$mean = rep(0,p)
    info$multiplier = diag(multiplier)

    result = list()
    result$pX = tmpdata%*%diag(multiplier)
    result$info = info
    return(result)
  } else if (type=="cscale"){ # add 2 : "center" + "scale"
    tmpdata = t(matinput)
    p       = ncol(tmpdata)

    infomean = as.double(colMeans(tmpdata))
    infomultiplier = rep(0,p)
    for (i in 1:p){
      infomultiplier[i] = 1/stats::sd(as.vector(tmpdata[,i]))
    }
    info = list()
    info$type = "cscale"
    info$mean = infomean
    info$multiplier = diag(infomultiplier)

    outdata = array(0,c(nrow(tmpdata),p)) # mean substraction
    for (i in 1:nrow(tmpdata)){
      outdata[i,] = as.vector(tmpdata[i,])-infomean
    }

    result = list()
    result$pX   = (outdata%*%diag(infomultiplier))
    result$info = info
    return(result)
  } else {
    # original methods
    # 'center','decorrelate', or 'whiten'
    matoutput = tryCatch(
      {
        if (type=="center"){
          aux_preprocess(matinput,as.integer(1))
        } else if (type=="decorrelate"){
          aux_preprocess(matinput,as.integer(2))
        } else {
          aux_preprocess(matinput,as.integer(3))
        }
      }, error=function(cond){
        return(NA)
      }, warning=function(cond){
        return(NA)
      }
    )
    if (length(matoutput)==1){
      if (is.na(matoutput)){
        result = NA
        return(result)
      }
    } else {
      # now we have 2 lists
      info   = list()
      info$type = matoutput$type
      info$mean = matoutput$mean
      info$multiplier = matoutput$multiplier

      result = list()
      result$pX = t(matoutput$output)
      result$info = info
      return(result)
    }
  }
}
#' @keywords internal
#' @noRd
aux.preprocess.hidden <- function(data,type=c("null","center","scale","cscale","decorrelate","whiten"),algtype=c("linear","nonlinear")){
  ## minimal preprocessing
  pptype    = match.arg(type)
  ppalgtype = match.arg(algtype)

  ## null is mine !
  n = nrow(data)
  p = ncol(data)
  if (type=="null"){
    info = list()
    info$type = "null"
    info$mean = rep(0,p)
    info$multiplier = 1

    result      = list()
    result$pX   = data
    result$info = info
    result$info$algtype = ppalgtype
    return(result)
  } else {
    result = aux.preprocess(data,type=pptype)
    result$info$algtype = ppalgtype
    return(result)
  }
}
