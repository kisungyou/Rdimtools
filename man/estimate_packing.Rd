% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_packing.R
\name{est.packing}
\alias{est.packing}
\title{Intrinsic Dimension Estimation using Packing Numbers}
\usage{
est.packing(X, eps = 0.01)
}
\arguments{
\item{X}{an \eqn{(n\times p)} matrix or data frame whose rows are observations.}

\item{eps}{small positive number for stopping threshold.}
}
\value{
a named list containing containing \describe{
\item{estdim}{estimated intrinsic dimension.}
}
}
\description{
Instead of covering numbers which are expensive to compute in many fractal-based methods,
\code{est.packing} exploits packing numbers as a proxy to describe spatial density. Since
it involves random permutation of the dataset at each iteration, every run might have
different results.
}
\examples{
\donttest{
## create 'swiss' roll dataset
X = aux.gensamples(dname="swiss")

## try different eps values
out1 = est.packing(X, eps=0.1)
out2 = est.packing(X, eps=0.01)
out3 = est.packing(X, eps=0.001)

## print the results
line1 = paste0("* est.packing : eps=0.1   gives ",round(out1$estdim,2))
line2 = paste0("* est.packing : eps=0.01  gives ",round(out2$estdim,2))
line3 = paste0("* est.packing : eps=0.001 gives ",round(out3$estdim,2))
cat(paste0(line1,"\n",line2,"\n",line3))
}

}
\references{
\insertRef{kegl_intrinsic_2002}{Rdimtools}
}
\author{
Kisung You
}
