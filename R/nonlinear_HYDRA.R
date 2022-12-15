#' Hyperbolic Distance Recovery and Approximation
#'
#' Hyperbolic Distance Recovery and Approximation, also known as \code{hydra} in short,
#' implements embedding of distance-based data into hyperbolic space represented as the Poincare disk,
#' which is interior of a hypersphere.
#'
#' @param X an \eqn{(n\times p)} matrix or data frame whose rows are observations and columns represent independent variables.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param ... extra parameters including \describe{
#' \item{kappa}{embedding curvature, which is a nonnegative number (default: 1).}
#' \item{iso.adjust}{perform isotropic adjustment. If \code{ndim=2}, default is \code{FALSE}. Otherwise, \code{TRUE} is used as default.}
#' }
#'
#' @return a named \code{Rdimtools} S3 object containing
#' \describe{
#' \item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations in the Poincare disk.}
#' \item{algorithm}{name of the algorithm.}
#' }
#'
#' @examples
#' \donttest{
#' ## load iris data
#' data(iris)
#' X     = as.matrix(iris[,1:4])
#' lab   = as.factor(iris[,5])
#'
#' ## multiple runs with varying curvatures
#' embed1 <- do.hydra(X, kappa=0.1)
#' embed2 <- do.hydra(X, kappa=1)
#' embed3 <- do.hydra(X, kappa=10)
#'
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(embed1$Y , col=lab, pch=19, main="kappa=0.1")
#' plot(embed2$Y , col=lab, pch=19, main="kappa=1")
#' plot(embed3$Y , col=lab, pch=19, main="kappa=10")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{keller-ressel_2020_HydraMethodStrainminimizing}{Rdimtools}
#'
#' @rdname nonlinear_HYDRA
#' @concept nonlinear_methods
#' @export
do.hydra <- function(X, ndim=2, ...){
  #------------------------------------------------------------------------
  # Preprocessing
  aux.typecheck(X)
  n = base::nrow(X)
  p = base::ncol(X)
  myndim = as.integer(ndim)
  if (!check_ndim(myndim,p)){
    stop("* do.hydra : 'ndim' is a positive integer in [1,#(covariates)].")
  }

  # Extra parameters
  params  = list(...)
  pnames  = names(params)
  if ("kappa" %in% pnames){
    kappa = max(sqrt(.Machine$double.eps), as.double(params$kappa))
  } else {
    kappa = 1
  }
  if ("iso.adjust" %in% pnames){
    iso_adjust = as.logical(params$iso.adjust)
  } else {
    if (myndim==2){
      iso_adjust = FALSE
    } else {
      iso_adjust = TRUE
    }
  }

  #------------------------------------------------------------------------
  # COMPUTATION
  # compute the pairwise distance
  DX = stats::dist(X)
  run_hydra   = my_hidden_hydra(DX, ndim=myndim, kappa=kappa, iso.adjust=iso_adjust)
  hydra_embed = array(0,c(n, myndim))
  for (i in 1:n){
    hydra_embed[i,] = run_hydra$radial[i]*as.vector(run_hydra$directional[i,])
  }

  #------------------------------------------------------------------------
  # Return
  result = list()
  result$Y = hydra_embed
  result$algorithm = "nonlinear:HYDRA"
  return(structure(result, class="Rdimtools"))
}



# auxiliary ---------------------------------------------------------------
#' @keywords internal
#' @noRd
my_hidden_hydra <- function(distobj, ndim=2, kappa=1, iso.adjust=TRUE){
  # preprocess
  D     = as.matrix(distobj)
  n     = base::nrow(D)
  dim   = round(ndim)
  kappa = base::max(base::sqrt(.Machine$double.eps), as.double(kappa))
  A     = base::cosh(base::sqrt(kappa)*D)

  # eigen-decomposition : 100 dimensions
  if (n < 100){
    # regular 'base::eigen'
    eigA = base::eigen(A, TRUE)

    # top vector
    x0 = sqrt(eigA$values[1])*as.vector(eigA$vectors[,1])
    if (x0[1] < 0){
      x0 = -x0
    }

    # others
    bot_vals = eigA$values[(n-dim+1):n]
    bot_vecs = eigA$vectors[,(n-dim+1):n]
  } else {
    # or use 'RSpectra::eigs_sym'
    eig_top = RSpectra::eigs_sym(A, 1,    which="LA")
    eig_bot = RSpectra::eigs_sym(A, ndim, which="SA")

    # top vector
    x0 = sqrt(eig_top$values)*as.vector(eig_top$vectors)
    if (x0[1] < 0){
      x0 = -x0
    }
    # others
    bot_vecs = eig_bot$vectors
    bot_vals = eig_bot$values
  }

  # component : radial
  x_min  = base::min(x0)
  radial = sqrt((x0-x_min)/(x0+x_min))

  # component : directional
  if (iso.adjust){
    X_last      = bot_vecs%*%diag(sqrt(pmax(-bot_vals,0)))
    sqnorms     = base::apply(X_last, 1, function(x) 1/sqrt(sum(x^2)))
    directional = base::diag(sqnorms)%*%X_last
  } else {
    sqnorms     = base::apply(bot_vecs, 1, function(x) 1/sqrt(sum(x^2)))
    directional = base::diag(sqnorms)%*%bot_vecs
  }

  # return the output
  return(list(radial=radial, directional=directional))
}
