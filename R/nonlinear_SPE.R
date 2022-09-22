#' Stochastic Proximity Embedding
#'
#' One of drawbacks for Multidimensional Scaling or Sammon mapping is that
#' they have quadratic computational complexity with respect to the number of data.
#' Stochastic Proximity Embedding (SPE) adopts stochastic update rule in that
#' its computational speed is much improved. It performs \code{C} number of cycles,
#' where for each cycle, it randomly selects two data points and updates their
#' locations correspondingly \code{S} times. After each cycle, learning parameter \eqn{\lambda}
#' is multiplied by \code{drate}, becoming smaller in magnitude.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations
#' and columns represent independent variables.
#' @param ndim an integer-valued target dimension.
#' @param proximity a function for constructing proximity matrix from original data dimension.
#' @param C the number of cycles to be run; after each cycle, learning parameter
#' @param S the number of updates for each cycle.
#' @param lambda initial learning parameter.
#' @param drate multiplier for \code{lambda} at each cycle; should be a positive real number in \eqn{(0,1).}
#'
#' @return a named list containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{trfinfo}{a list containing information for out-of-sample prediction.}
#' }
#'
#' @examples
#' \donttest{
#' ## load iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' label = as.factor(iris$Species)
#'
#' ## compare with mds using 2 distance metrics
#' outM <- do.mds(X, ndim=2)
#' out1 <- do.spe(X, ndim=2)
#' out2 <- do.spe(X, ndim=2, proximity=function(x){dist(x, method="manhattan")})
#'
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(outM$Y, pch=19, col=label, main="MDS")
#' plot(out1$Y, pch=19, col=label, main="SPE with L2 norm")
#' plot(out2$Y, pch=19, col=label, main="SPE with L1 norm")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{agrafiotis_stochastic_2003}{Rdimtools}
#'
#' @author Kisung You
#' @rdname nonlinear_SPE
#' @concept nonlinear_methods
#' @export
do.spe <- function(X, ndim=2, proximity=function(x){dist(x, method="euclidean")},
                   C = 50, S = 50, lambda = 1, drate=0.9){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  # 1. data X
  aux.typecheck(X)
  n = nrow(X)
  p = ncol(X)
  # 2. ndim
  if ((!is.numeric(ndim))||(ndim<1)||(ndim>=ncol(X))||is.infinite(ndim)||is.na(ndim)){
    stop("* do.spe : 'ndim' is a positive integer in [1,#(covariates)).")
  }
  # 3. proximity
  if (!is.function(proximity)){
    stop("* do.spe : 'proximity' should be a function for computing proximity measure.")
  }
  # 4. C
  C = as.integer(C)
  if (!check_NumMM(C,2,Inf,compact=TRUE)){stop("* do.spe : 'C' should be a positive integer number.")}
  # 5. S
  S = as.integer(S)
  if (!check_NumMM(S,2,Inf,compact=TRUE)){stop("* do.spe : 'S' should be a positive integer number.")}
  # 6. lambda
  lambda = as.double(lambda)
  if (!check_NumMM(lambda,0,Inf, compact=FALSE)){stop("* do.spe : 'lambda' should be a positive real number.")}
  # 7. drate
  drate = as.double(drate)
  if (!check_NumMM(drate,0,1,compact=FALSE)){stop("* do.spe : 'drate' should be in (0,1).")}


  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  # 1. Proximity Matrix R : don't forget to make it into a matrix.
  R = proximity(X)
  if (inherits(R, "dist")){
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
  Y = method_spe(R, iX, C, S, lambda, drate, matselector)


  #------------------------------------------------------------------------
  ## RETURN OUTPUT
  trfinfo = list()
  trfinfo$type = "null"
  trfinfo$algtype = "nonlinear"
  trfinfo$mean = rep(0,p)
  trfinfo$multiplier = 1

  result = list()
  result$Y = Y
  result$trfinfo = trfinfo
  return(result)
}
