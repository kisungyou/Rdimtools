% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linear_SPP.R
\name{do.spp}
\alias{do.spp}
\title{Sparsity Preserving Projection}
\usage{
do.spp(
  X,
  ndim = 2,
  preprocess = c("center", "scale", "cscale", "decorrelate", "whiten"),
  reltol = 1e-04
)
}
\arguments{
\item{X}{an \eqn{(n\times p)} matrix or data frame whose rows are observations}

\item{ndim}{an integer-valued target dimension.}

\item{preprocess}{an additional option for preprocessing the data.
Default is "center". See also \code{\link{aux.preprocess}} for more details.}

\item{reltol}{tolerance level for stable computation of sparse reconstruction weights.}
}
\value{
a named list containing
\describe{
\item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
\item{trfinfo}{a list containing information for out-of-sample prediction.}
\item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
}
}
\description{
Sparsity Preserving Projection (SPP) is an unsupervised linear dimension reduction technique.
It aims to preserve high-dimensional structure in a sparse manner to find projections
that keeps such sparsely-connected pattern in the low-dimensional space. Note that
we used \pkg{CVXR} for convenient computation, which may lead to slower execution
once used for large dataset.
}
\examples{
\dontrun{
## load iris data
data(iris)
set.seed(100)
subid = sample(1:150,50)
X     = as.matrix(iris[subid,1:4])
label = as.factor(iris[subid,5])

## test different tolerance levels
out1 <- do.spp(X,ndim=2,reltol=0.001)
out2 <- do.spp(X,ndim=2,reltol=0.01)
out3 <- do.spp(X,ndim=2,reltol=0.1)

# visualize
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3))
plot(out1$Y, pch=19, col=label, main="SPP::reltol=.001")
plot(out2$Y, pch=19, col=label, main="SPP::reltol=.01")
plot(out3$Y, pch=19, col=label, main="SPP::reltol=.1")
par(opar)
}
}
\references{
\insertRef{qiao_sparsity_2010}{Rdimtools}
}
\author{
Kisung You
}
\concept{linear_methods}
