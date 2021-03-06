% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_MIFS.R
\name{do.mifs}
\alias{do.mifs}
\title{Mutual Information for Selecting Features}
\usage{
do.mifs(
  X,
  label,
  ndim = 2,
  beta = 0.75,
  discretize = c("default", "histogram"),
  preprocess = c("null", "center", "scale", "cscale", "whiten", "decorrelate")
)
}
\arguments{
\item{X}{an \eqn{(n\times p)} matrix or data frame whose rows are observations
and columns represent independent variables.}

\item{label}{a length-\eqn{n} vector of class labels.}

\item{ndim}{an integer-valued target dimension.}

\item{beta}{penalty for relative importance of mutual information between the candidate and already-chosen features in iterations. Author proposes to use a value in \eqn{(0.5,1)}.}

\item{discretize}{the method for each variable to be discretized. The paper proposes \code{"default"} method to use 10 bins while \code{"histogram"} uses automatic discretization via Sturges' method.}

\item{preprocess}{an additional option for preprocessing the data. Default is "null". See also \code{\link{aux.preprocess}} for more details.}
}
\value{
a named list containing
\describe{
\item{Y}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
\item{featidx}{a length-\eqn{ndim} vector of indices with highest scores.}
\item{trfinfo}{a list containing information for out-of-sample prediction.}
\item{projection}{a \eqn{(p\times ndim)} whose columns are basis for projection.}
}
}
\description{
MIFS is a supervised feature selection that iteratively increases the subset of variables by choosing maximally informative feature based on the mutual information.
}
\examples{
\donttest{
## use iris data
## it is known that feature 3 and 4 are more important.
data(iris)
iris.dat = as.matrix(iris[,1:4])
iris.lab = as.factor(iris[,5])

## try different beta values
out1 = do.mifs(iris.dat, iris.lab, beta=0)
out2 = do.mifs(iris.dat, iris.lab, beta=0.5)
out3 = do.mifs(iris.dat, iris.lab, beta=1)

## visualize
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3))
plot(out1$Y, pch=19, col=iris.lab, main="beta=0")
plot(out2$Y, pch=19, col=iris.lab, main="beta=0.5")
plot(out3$Y, pch=19, col=iris.lab, main="beta=1")
par(opar)
}


}
\references{
\insertRef{battiti_using_1994}{Rdimtools}
}
\author{
Kisung You
}
\concept{feature_methods}
