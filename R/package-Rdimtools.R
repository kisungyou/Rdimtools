#' Dimension Reduction and Estimation Methods
#'
#' \pkg{Rdimtools} is an R implementation of a number of dimension reduction and estimation methods
#' implemented using \pkg{RcppArmadillo} for efficient computations. Please see the section below for
#' the complete composition of this package and what we can provide in a unifying interface across many
#' methods.
#'
#' @section Composition of the package:
#' The package consists of three families of functions whose names start with \code{do.}, \code{est.}, and \code{aux.}
#' for performing dimension reduction/manifold learning, estimating intrinsic dimension, and some efficient
#' implementations of other useful methods respectively.
#'
#' \bold{\code{do.}} functions are for dimension reduction (or, \emph{manifold learning}) methods.
#' A simple taxonomy of the methods would be to categorize based on the linearity of
#' embedding mappings. In the table below, TYPE represents whether it is \emph{supervised}(S),
#' \emph{semisupervised}(SS), or \emph{unsupervised}(U).
#'
#' For \emph{linear} methods, we have
#' \tabular{lcl}{
#' FUNCTION \tab TYPE \tab ALGORITHM \cr
#' \code{\link{do.anmm}}\tab S \tab Average Neighborhood Margin Maximization \cr
#' \code{\link{do.bpca}} \tab U \tab Bayesian Principal Component Analysis \cr
#' \code{\link{do.cca}} \tab S \tab Canonical Correlation Analysis \cr
#' \code{\link{do.fa}} \tab U \tab (Exploratory) Factor Analysis  \cr
#' \code{\link{do.ica}} \tab U \tab Independent Component Analysis \cr
#' \code{\link{do.isoproj}} \tab U \tab Isometric projection \cr
#' \code{\link{do.lda}} \tab S \tab Linear Discriminant Analysis \cr
#' \code{\link{do.lmds}} \tab U \tab Landmark Multidimensional Scaling \cr
#' \code{\link{do.lpp}} \tab U \tab Locality Preserving Embedding \cr
#' \code{\link{do.mds}} \tab U \tab (Metric) Multidimensional Scaling \cr
#' \code{\link{do.npe}} \tab U \tab Neighborhood Preserving Embedding \cr
#' \code{\link{do.olpp}} \tab U \tab Orthogonal Locality Preserving Embedding \cr
#' \code{\link{do.opls}} \tab S \tab Orthogonal Partial Least Squares \cr
#' \code{\link{do.pca}} \tab U \tab Principal Component Analysis \cr
#' \code{\link{do.pls}} \tab S \tab Partial Least Squares \cr
#' \code{\link{do.ppca}} \tab U \tab Probabilistic Principal Component Analysis \cr
#' \code{\link{do.rndproj}} \tab U \tab Random Projection \cr
#' \code{\link{do.spca}} \tab U \tab Sparse Principal Component Analysis
#' }
#'
#' Also, we have \emph{nonlinear} methods implemented
#' \tabular{lcl}{
#' FUNCTION \tab TYPE \tab ALGORITHM \cr
#' \code{\link{do.cisomap}} \tab U \tab Conformal Isometric Feature Mapping \cr
#' \code{\link{do.dm}} \tab U \tab Diffusion Maps \cr
#' \code{\link{do.isomap}} \tab U \tab Isometric Feature Mapping \cr
#' \code{\link{do.ispe}} \tab U \tab Isometric Stochastic Proximity Embedding \cr
#' \code{\link{do.keca}} \tab U \tab Kernel Entropy Component Analysis \cr
#' \code{\link{do.kpca}} \tab U \tab Kernel Principal Component Analysis \cr
#' \code{\link{do.lapeig}} \tab U \tab Laplacian Eigenmaps \cr
#' \code{\link{do.lisomap}} \tab U \tab Landmark Isometric Feature Mapping \cr
#' \code{\link{do.lle}} \tab U \tab Locally Linear Embedding \cr
#' \code{\link{do.ltsa}} \tab U \tab Local Tangent Space Alignment \cr
#' \code{\link{do.mve}} \tab U \tab Minimum Volume Embedding \cr
#' \code{\link{do.mvu}} \tab U \tab Maximum Variance Unfolding / Semidefinite Embedding \cr
#' \code{\link{do.plp}} \tab U \tab Piecewise Laplacian Projection \cr
#' \code{\link{do.ree}} \tab U \tab Robust Euclidean Embedding \cr
#' \code{\link{do.rpca}}\tab U \tab Robust Principal Component Analysis \cr
#' \code{\link{do.sammon}} \tab U \tab Sammon Mapping \cr
#' \code{\link{do.sne}} \tab U \tab Stochastic Neighbor Embedding \cr
#' \code{\link{do.spe}} \tab U \tab Stochastic Proximity Embedding \cr
#' \code{\link{do.tsne}} \tab U \tab t-distributed Stochastic Neighbor Embedding
#' }
#'
#' Secondly, \bold{\code{est.}} family of functions are for intrinsic dimension estimation methods, including
#' \tabular{ll}{
#' FUNCTION \tab ALGORITHM \cr
#' \code{\link{est.boxcount}} \tab Box-Counting Dimension \cr
#' \code{\link{est.correlation}} \tab Correlation Dimension
#' }
#'
#' Finally, there are some auxiliary functions (\bold{\code{aux.}} family),
#' \tabular{ll}{
#' FUNCTION \tab DESCRIPTION \cr
#' \code{\link{aux.gensamples}} \tab generate samples from predefined shapes \cr
#' \code{\link{aux.graphnbd}} \tab make a neighborhood graph given certain criteria \cr
#' \code{\link{aux.kernelcov}} \tab compute a centered gram matrix with 20 kernels supported \cr
#' \code{\link{aux.preprocess}} \tab perform preprocessing of centering, decorrelating, or whitening \cr
#' \code{\link{aux.shortestpath}} \tab Floyd-Warshall algorithm (it's \emph{Fast}!) \cr
#' \code{\link{aux.pkgstat}} \tab report the number of functions available for each category
#' }
#'
#'
#' @docType package
#' @name Rdimtools
#' @aliases Rdimtools-package
#' @import Rcsdp
#' @import Rdpack
#' @import RSpectra
#' @import CVXR
#' @importFrom Rlinsolve lsolve.bicgstab
#' @importFrom Rtsne Rtsne
#' @importFrom stats dist cov rnorm runif
#' @importFrom graphics par image plot
#' @importFrom Rcpp evalCpp
#' @useDynLib Rdimtools
NULL

# NOTES
# 1. reference travis builder : https://github.com/mlr-org/mlr/blob/master/.travis.yml
