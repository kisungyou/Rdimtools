#' Isometric Stochastic Proximity Embedding
#'
#' The isometric SPE (ISPE) adopts the idea of approximating geodesic distance on embedded manifold
#' when two data points are close enough. It introduces the concept of \code{cutoff} where the learning process
#' is only applied to the pair of data points whose original proximity is small enough to be considered as
#' mutually local whose distance should be close to geodesic distance.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param proximity a function for constructing proximity matrix from original data dimension.
#' @param C the number of cycles to be run; after each cycle, learning parameter
#' @param S the number of updates for each cycle.
#' @param lambda initial learning parameter.
#' @param drate multiplier for \code{lambda} at each cycle; should be a positive real number in \eqn{(0,1).}
#' @param cutoff cutoff threshold value.
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' ## generate sample data
#' X = aux.gensamples()
#'
#' ## compare with original SPE
#' outSPE <- do.spe(X, ndim=2)
#' out1 <- do.ispe(X, ndim=2, cutoff=0.5)
#' out2 <- do.ispe(X, ndim=2, cutoff=5)
#' out3 <- do.ispe(X, ndim=2, cutoff=50)
#'
#' ## Visualize
#' par(mfrow=c(2,2))
#' plot(outSPE$Y[,1], outSPE$Y[,2], main="SPE")
#' plot(out1$Y[,1], out1$Y[,2], main="ISPE with cutoff=0.5")
#' plot(out2$Y[,1], out2$Y[,2], main="ISPE with cutoff=5")
#' plot(out3$Y[,1], out3$Y[,2], main="ISPE with cutoff=50")
#'
#' @references
#' \insertRef{agrafiotis_self-organizing_2002}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_ISPE
#' @export
do.ispe <- function(X, ndim=2, proximity=function(x){dist(x, method="euclidean")},
                   C = 50, S = 50, lambda = 1, drate=0.9, cutoff = 1.0){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  # 1. data X
  aux.typecheck(X)
  # 2. ndim
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>=ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.ispe : 'ndim' is a positive integer in [1,#(covariates)).")
  }
  # 3. proximity
  if (!is.function(proximity)){
    stop("* do.ispe : 'proximity' should be a function for computing proximity measure.")
  }
  # 4. C
  C = as.integer(C)
  if (!check_NumMM(C,2,Inf,compact=TRUE)){stop("* do.ispe : 'C' should be a positive integer number.")}
  # 5. S
  S = as.integer(S)
  if (!check_NumMM(S,2,Inf,compact=TRUE)){stop("* do.ispe : 'S' should be a positive integer number.")}
  # 6. lambda
  lambda = as.double(lambda)
  if (!check_NumMM(lambda,0,Inf, compact=FALSE)){stop("* do.ispe : 'lambda' should be a positive real number.")}
  # 7. drate
  drate = as.double(drate)
  if (!check_NumMM(drate,0,1,compact=FALSE)){stop("* do.ispe : 'drate' should be in (0,1).")}
  # 8. cutoff
  cutoff = as.double(cutoff)
  if (!check_NumMM(cutoff,0,Inf, compact=FALSE)){stop("* do.ispe : 'cutoff' should be a positive real number.")}


  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  # 1. Proximity Matrix R : don't forget to make it into a matrix.
  R = proximity(X)
  if ((inherits(R, "dist"))||(class(R)=="dist")){
    R = as.matrix(R)
  }
  # 2. initial coordinates
  iX = do.pca(X, ndim=ndim)$Y
  # 3. index vector
  vecselector = seq(from=0,to=(nrow(X)-1),by=1)
  matselector = array(0,c(C*S,2))
  for (i in 1:(C*S)){
    matselector[i,] = sample(vecselector,2)
  }
  # 4. pass onto Rcpp
  Y = method_ispe(R, iX, C, S, lambda, drate, matselector, cutoff)


  #------------------------------------------------------------------------
  ## RETURN OUTPUT
  trfinfo = list()
  trfinfo$type = "null"
  trfinfo$algtype = "nonlinear"

  result = list()
  result$Y = Y
  result$trfinfo = trfinfo
  return(result)
}
