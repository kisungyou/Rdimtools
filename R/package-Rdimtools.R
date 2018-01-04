#' Dimension Reduction and Estimation Methods
#'
#' \pkg{Rdimtools} is an R suite of a number of dimension reduction and estimation methods
#' implemented using \pkg{RcppArmadillo} for efficient computations. Please see the section below for
#' the complete composition of this package and what we can provide in a unifying interface across many
#' methods.
#'
#' @section Composition of the package:
#' The package consists of following families of functions whose names start with \code{do.}, \code{est.}, and \code{aux.}
#' for performing dimension reduction/manifold learning, estimating intrinsic dimension, and some efficient
#' implementations of other useful methods respectively.
#'
#' @section (1) \code{do.} family for dimension reduction algorithms:
#' \bold{\code{do.}} functions are for dimension reduction (or, \emph{manifold learning}) methods.
#' A simple taxonomy of the methods would be to categorize based on the linearity of
#' embedding mappings. In the table below, TYPE represents whether it is \emph{supervised} (S),
#' \emph{semisupervised} (SS), or \emph{unsupervised} (U).
#'
#' For \emph{linear} methods, we have
#' \tabular{lcl}{
#' FUNCTION \tab TYPE \tab ALGORITHM \cr
#' \code{\link{do.anmm}}\tab S \tab Average Neighborhood Margin Maximization \cr
#' \code{\link{do.bpca}} \tab U \tab Bayesian Principal Component Analysis \cr
#' \code{\link{do.cca}} \tab S \tab Canonical Correlation Analysis \cr
#' \code{\link{do.dspp}} \tab S \tab Discriminative Sparsity Preserving Projection \cr
#' \code{\link{do.eslpp}} \tab S \tab Extended Supervised Locality Preserving Projection \cr
#' \code{\link{do.extlpp}} \tab U \tab Extended Locality Preserving Projection \cr
#' \code{\link{do.fa}} \tab U \tab (Exploratory) Factor Analysis  \cr
#' \code{\link{do.ica}} \tab U \tab Independent Component Analysis \cr
#' \code{\link{do.isoproj}} \tab U \tab Isometric Projection \cr
#' \code{\link{do.kmvp}} \tab S \tab Kernel-Weighted Maximum Variance Projection \cr
#' \code{\link{do.kudp}} \tab U \tab Kernel-Weighted Unsupervised Discriminant Projection \cr
#' \code{\link{do.lda}} \tab S \tab Linear Discriminant Analysis \cr
#' \code{\link{do.lde}} \tab S \tab Local Discriminant Embedding \cr
#' \code{\link{do.lea}} \tab U \tab Locally Linear Embedded Eigenspace Analysis \cr
#' \code{\link{do.lfda}} \tab S \tab Local Fisher Discriminant Analysis \cr
#' \code{\link{do.llp}} \tab U \tab Local Learning Projections \cr
#' \code{\link{do.lmds}} \tab U \tab Landmark Multidimensional Scaling \cr
#' \code{\link{do.lpp}} \tab U \tab Locality Preserving Projection \cr
#' \code{\link{do.lspp}}\tab S \tab Local Similarity Preserving Projection \cr
#' \code{\link{do.mds}} \tab U \tab (Metric) Multidimensional Scaling \cr
#' \code{\link{do.mmc}} \tab S \tab Maximum Margin Criterion \cr
#' \code{\link{do.modp}} \tab S \tab Modified Orthogonal Discriminant Projection \cr
#' \code{\link{do.mvp}} \tab S \tab Maximum Variance Projection \cr
#' \code{\link{do.npe}} \tab U \tab Neighborhood Preserving Embedding \cr
#' \code{\link{do.odp}} \tab S \tab Orthogonal Discriminant Projection \cr
#' \code{\link{do.olpp}} \tab U \tab Orthogonal Locality Preserving Projection \cr
#' \code{\link{do.opls}} \tab S \tab Orthogonal Partial Least Squares \cr
#' \code{\link{do.pca}} \tab U \tab Principal Component Analysis \cr
#' \code{\link{do.pls}} \tab S \tab Partial Least Squares \cr
#' \code{\link{do.ppca}} \tab U \tab Probabilistic Principal Component Analysis \cr
#' \code{\link{do.rda}} \tab S \tab Regularized Discriminant Analysis \cr
#' \code{\link{do.rndproj}} \tab U \tab Random Projection \cr
#' \code{\link{do.sda}} \tab SS \tab Semi-Supervised Discriminant Analysis \cr
#' \code{\link{do.sdlpp}} \tab U \tab Sample-Dependent Locality Preserving Projection \cr
#' \code{\link{do.slpp}} \tab S \tab Supervised Locality Preserving Projection \cr
#' \code{\link{do.spca}} \tab U \tab Sparse Principal Component Analysis \cr
#' \code{\link{do.spp}} \tab U \tab Sparsity Preserving Projection \cr
#' \code{\link{do.udp}} \tab U \tab Unsupervised Discriminant Projection
#' }
#'
#' Also, we have \emph{nonlinear} methods implemented
#' \tabular{lcl}{
#' FUNCTION \tab TYPE \tab ALGORITHM \cr
#' \code{\link{do.cisomap}} \tab U \tab Conformal Isometric Feature Mapping \cr
#' \code{\link{do.crca}} \tab U \tab Curvilinear Component Analysis \cr
#' \code{\link{do.crda}} \tab U \tab Curvilinear Distance Analysis \cr
#' \code{\link{do.dm}} \tab U \tab Diffusion Maps \cr
#' \code{\link{do.isomap}} \tab U \tab Isometric Feature Mapping \cr
#' \code{\link{do.ispe}} \tab U \tab Isometric Stochastic Proximity Embedding \cr
#' \code{\link{do.keca}} \tab U \tab Kernel Entropy Component Analysis \cr
#' \code{\link{do.klde}} \tab S \tab Kernel Local Discriminant Embedding \cr
#' \code{\link{do.klfda}} \tab S \tab Kernel Local Fisher Discriminant Analysis \cr
#' \code{\link{do.kmmc}} \tab S \tab Kernel Maximium Margin Criterion \cr
#' \code{\link{do.kpca}} \tab U \tab Kernel Principal Component Analysis \cr
#' \code{\link{do.ksda}} \tab SS \tab Kernel Semi-Supervised Discriminant Analysis \cr
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
#' @section (2) \code{est.} family for intrinsic dimension estimation algorithms:
#' \bold{\code{est.}} family of functions include,
#' \tabular{ll}{
#' FUNCTION \tab ALGORITHM \cr
#' \code{\link{est.boxcount}} \tab Box-Counting Dimension \cr
#' \code{\link{est.correlation}} \tab Correlation Dimension
#' }
#'
#' @section (3) \code{oos.} family for out-of-sample predictions:
#' Only except for \code{oos.linear} where we have explicit mappings from linear dimension reduction methods, all methods
#' listed in this section are for out-of-sample prediction on nonlinear mappings.
#' \tabular{lcl}{
#' FUNCTION \tab LINEARITY \tab ALGORITHM \cr
#' \code{\link{oos.linear}} \tab Linear \tab Projection using pre-computed mappings.
#' }
#'
#'
#' @section (4) \code{aux.} functions:
#' Some auxiliary functions (\bold{\code{aux.}}) are also provided,
#' \tabular{ll}{
#' FUNCTION \tab DESCRIPTION \cr
#' \code{\link{aux.gensamples}} \tab generates samples from predefined shapes \cr
#' \code{\link{aux.graphnbd}} \tab builds a neighborhood graph given certain criteria \cr
#' \code{\link{aux.kernelcov}} \tab computes a centered gram matrix with 20 kernels supported \cr
#' \code{\link{aux.preprocess}} \tab performs preprocessing of centering, decorrelating, or whitening \cr
#' \code{\link{aux.shortestpath}} \tab Floyd-Warshall algorithm (it's \emph{Fast}!) \cr
#' \code{\link{aux.pkgstat}} \tab reports the number of functions available for each category
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
#' @import geigen
#' @importFrom Matrix rankMatrix
#' @importFrom Rlinsolve lsolve.bicgstab
#' @importFrom Rtsne Rtsne
#' @importFrom stats dist cov rnorm runif kmeans cor
#' @importFrom graphics par image plot
#' @importFrom Rcpp evalCpp
#' @useDynLib Rdimtools
NULL

# NOTES
# 1. reference travis builder : https://github.com/mlr-org/mlr/blob/master/.travis.yml
