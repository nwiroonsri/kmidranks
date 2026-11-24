\name{Acc_label}
\alias{Acc_label}
\title{Clustering Accuracy After Optimal Label Matching}
\description{Calculate the largest fraction of observations correctly
classified under a one-to-one matching of cluster labels to reference labels.}
\usage{
Acc_label(kmc, y)
}
\arguments{
\item{kmc}{Atomic vector of predicted cluster labels, such as the
\code{cluster} component of a \code{\link{K_midranks}} fit.}
\item{y}{Atomic vector of reference labels for the same observations in the
same order. Both vectors must be nonempty, equally long, and free of missing values.}
}
\details{
The counts of predicted labels by reference labels form a contingency table.
A linear assignment maximizes the sum of matched counts. If there are more
predicted clusters than reference classes, extra clusters can match zero-count
dummy columns. If there are fewer, some reference classes remain unmatched.
Each unmatched observation contributes zero to the number correctly classified.
}
\value{A single numeric accuracy between zero and one.}
\examples{
data("cauchy_k3", package = "Kmidranks")
x <- subset(cauchy_k3, select = -label)
set.seed(42)
fit <- K_midranks(x, k = 3, nstart = 3)
table(cluster = fit$cluster, truth = cauchy_k3$label)
Acc_label(fit$cluster, cauchy_k3$label)
}
\seealso{\code{\link{K_midranks}}}
