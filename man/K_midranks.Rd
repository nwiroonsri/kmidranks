\name{K_midranks}
\alias{K_midranks}
\title{Cluster Observations Using Feature-Wise Ranks}
\description{Partition the rows of a numeric data set by alternating cluster
assignment and center updates in rank space.}
\usage{
K_midranks(x, k, dist = 1, nstart = 20, mr = 1000)
}
\arguments{
\item{x}{Numeric matrix or all-numeric data frame. Rows are observations and
columns are features. Requires at least two rows, one column, and finite entries.}
\item{k}{Integer number of clusters, from one through the number of rows. With one cluster, the center and objective are computed directly (zero iterations).}
\item{dist}{Rank-space dissimilarity: \code{1} for Manhattan distance,
\code{2} for squared Euclidean distance, or \code{Inf} for the maximum
absolute coordinate difference.}
\item{nstart}{Positive integer number of random starting partitions.}
\item{mr}{Positive integer maximum number of assignment/update iterations per start.}
}
\details{
Each rank equals one plus the number of observations with a strictly larger
value in that feature. For example, values 9, 7, 7, 2 receive ranks 1, 2, 2, 4.
A constant feature has rank 1 throughout and contributes zero distance;
such features trigger a warning.

Starting partitions randomly permute a balanced set of cluster labels.
Centers are coordinate-wise medians for \code{dist = 1}, means for
\code{dist = 2}, and maxima for \code{dist = Inf}. Equal assignment
distances favor the lowest cluster index. Empty clusters are repaired by
moving observations from clusters of size at least two. Iteration stops when
the repaired assignments are unchanged or the iteration limit is reached.
The first start attaining the smallest final total dissimilarity is retained.
Use \code{set.seed()} before fitting for reproducibility.
}
\value{A list describing the retained start:
\item{cluster}{Integer vector of cluster assignments in input row order.}
\item{center}{Matrix with one row per cluster and one column per feature,
in rank units.}
\item{size}{Integer vector of cluster sizes in cluster-index order.}
\item{within}{Vector of dissimilarities from each observation to its assigned center.}
\item{tot.within}{Sum of \code{within}.}
\item{iter}{Number of iterations performed in the retained start.}
\item{converged}{Logical indicating whether assignments stabilized.
A warning is issued if the retained start reached the limit without stabilizing.}
}
\examples{
data("cauchy_k3", package = "Kmidranks")
x <- subset(cauchy_k3, select = -label)
set.seed(42)
fit <- K_midranks(x, k = 3, nstart = 3)
table(cluster = fit$cluster, truth = cauchy_k3$label)
Acc_label(fit$cluster, cauchy_k3$label)
fit$size
fit$center
stopifnot(isTRUE(all.equal(sum(fit$within), fit$tot.within)))
}
\seealso{\code{\link{elbow_rk}}, \code{\link{Acc_label}}}
