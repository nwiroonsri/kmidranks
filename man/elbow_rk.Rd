\name{elbow_rk}
\alias{elbow_rk}
\title{Compare Cluster Counts Using Successive Objective Reductions}
\description{Fit K-midranks over a range of cluster counts and rank candidate
elbows by the ratio of the objective reductions on either side.}
\usage{
elbow_rk(x, dist = 1, kmax = 10, nstart = 5, mr = 1000,
         verbose = TRUE, plot = FALSE)
}
\arguments{
\item{x}{Numeric matrix or all-numeric data frame with finite entries, at
least one feature and at least \code{kmax + 1} rows.}
\item{dist}{Dissimilarity passed to \code{\link{K_midranks}}:
\code{1}, \code{2}, or \code{Inf}.}
\item{kmax}{Integer upper limit for candidate elbows, at least two. Fits
are calculated for cluster counts from one through \code{kmax + 1}.}
\item{nstart}{Positive integer number of random starts for each cluster count.}
\item{mr}{Positive integer iteration limit per start.}
\item{verbose}{Logical. Print progress as each cluster count is fitted?}
\item{plot}{Logical. If true, plot the total within-cluster dissimilarity
against the number of clusters. A dashed line marks \code{best.K} when a
finite elbow ratio is available.}
}
\details{
Let \eqn{W_k} denote the fitted total dissimilarity and
\eqn{d_k = W_k - W_{k+1}}. For each candidate from two through
\code{kmax}, the elbow ratio is \eqn{d_{k-1}/d_k}.
The one-cluster fit supplies the baseline for the first ratio,
\eqn{(W_1-W_2)/(W_2-W_3)}. A large ratio means the reduction on reaching that count was large relative
to the reduction from adding one more cluster. A ratio is recorded as
\code{NA} if either reduction is nonfinite, the numerator is negative,
or the denominator is nonpositive. Undefined ratios produce a warning.
Only finite ratios are eligible for selection; ties favor smaller counts.
A fit that fails to converge generates a warning and its objective is still used.
}
\value{A list with components:
\item{best.K, second.K}{Cluster counts with the largest and second-largest
eligible ratios. A missing selection is returned as \code{NA_integer_}.}
\item{ebv}{Matrix of cluster counts (one through \code{kmax}) and forward
objective reductions, with columns \code{K} and \code{slope}.}
\item{ebi}{Matrix of candidate counts and elbow ratios, with columns
\code{K} and \code{ratio}.}
\item{objective}{Matrix of all fitted counts and total dissimilarities,
with columns \code{K} and \code{within}.}
\item{converged}{Matrix with columns \code{K} and \code{converged};
the latter is numeric, one for convergence and zero otherwise.}
\item{elapsed}{Elapsed computation time in seconds.}
}
\examples{
data("cauchy_k3", package = "Kmidranks")
x <- subset(cauchy_k3, select = -label)
set.seed(42)
e <- elbow_rk(x, kmax = 5, nstart = 3, verbose = FALSE, plot = TRUE)
e$ebi
e$best.K
}
\seealso{\code{\link{K_midranks}}}
