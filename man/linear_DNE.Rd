% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linear_DNE.R
\name{do.dne}
\alias{do.dne}
\title{Discriminant Neighborhood Embedding}
\usage{
do.dne(
  X,
  label,
  ndim = 2,
  numk = max(ceiling(nrow(X)/10), 2),
  preprocess = c("center", "scale", "cscale", "decorrelate", "whiten")
)
}
\arguments{
\item{X}{an \eqn{(n\times p)} matrix or data frame whose rows are observations.}

\item{label}{a length-\eqn{n} vector of data class labels.}

\item{ndim}{an integer-valued target dimension.}

\item{numk}{the number of neighboring points for k-nn graph construction.}

\item{preprocess}{an additional option for preprocessing the data.
Default is "center". See also \code{\link{aux.preprocess}} for more details.}
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
Discriminant Neighborhood Embedding (DNE) is a supervised subspace learning method.
DNE tries to move multi-class data points in high-dimensional space in accordance with
local intra-class attraction and inter-class repulsion.
}
\examples{
## load iris data
data(iris)
set.seed(100)
subid = sample(1:150,50)
X     = as.matrix(iris[subid,1:4])
label = as.factor(iris[subid,5])

## try different numbers for neighborhood size
out1 = do.dne(X, label, numk=5)
out2 = do.dne(X, label, numk=10)
out3 = do.dne(X, label, numk=20)

## visualize
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3))
plot(out1$Y, main="DNE::nbd size=5",  col=label, pch=19)
plot(out2$Y, main="DNE::nbd size=10", col=label, pch=19)
plot(out3$Y, main="DNE::nbd size=20", col=label, pch=19)
par(opar)

}
\references{
\insertRef{zhang_discriminant_2006}{Rdimtools}
}
\author{
Kisung You
}
\concept{linear_methods}
