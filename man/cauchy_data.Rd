\name{cauchy_data}
\alias{cauchy_k2}
\alias{cauchy_k3}
\alias{cauchy_k4}
\alias{cauchy_k5}
\alias{cauchy_k6}
\alias{cauchy_k7}
\alias{cauchy_k8}
\alias{cauchy_k9}
\docType{data}
\title{Simulated Clusters with Group-Specific Cauchy Perturbations}
\description{Eight simulated datasets from the K-midranks clustering study,
one for each number of clusters from two through nine.}
\usage{
data(cauchy_k2)
data(cauchy_k3)
data(cauchy_k4)
data(cauchy_k5)
data(cauchy_k6)
data(cauchy_k7)
data(cauchy_k8)
data(cauchy_k9)
}
\format{For a dataset with \eqn{K} groups, the data frame has \eqn{1000 K}
rows and \eqn{2K + 1} columns. The first \eqn{2K} columns are numeric
features named \code{x1}, \code{x2}, and so on. The final column,
\code{label}, gives the true group as a numeric value from one through
\eqn{K}. Each group occupies a block of 1,000 consecutive rows.}
\details{All features begin as independent standard normal draws. For
observations in group \eqn{i}, independent Cauchy draws with location 2
and scale 0.3 are added to features \eqn{2i-1} and \eqn{2i}.
For example, group 2 receives these additions in \code{x3} and \code{x4}.
Exclude \code{label} when fitting a clustering algorithm.}
\examples{
data(cauchy_k3)
dim(cauchy_k3)
table(cauchy_k3$label)
}
\keyword{datasets}
