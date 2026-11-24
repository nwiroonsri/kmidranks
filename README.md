<p align="center">
  <img src="Kmidranks ranking clusters hex logo.png" width="200" />
</p>

# Kmidranks

K-midranks clusters observations using the ranks of their features. This is
useful for heavy-tailed measurements when two unusually large observations
should differ according to their ordering rather than the gap between their
raw values. Within each feature, the largest values receive rank 1.

## Installation

Install Kmidranks from GitHub in R:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("nwiroonsri/kmidranks")
library(Kmidranks)
```

### Installation from a downloaded source folder

Unzip the source folder. Install the dependency in R:

```r
install.packages("clue")
```

From a terminal in the directory containing the extracted Kmidranks folder:

```sh
R CMD INSTALL Kmidranks
```

## Fit and evaluate a partition

The package includes eight simulated datasets, `cauchy_k2` through
`cauchy_k9`. In `cauchy_k3`, each of the three groups has 1,000
observations and a different pair of features with added Cauchy noise.
Fit the six features and reserve the known labels for evaluation.

```r
library(Kmidranks)
data("cauchy_k3", package = "Kmidranks")
x <- subset(cauchy_k3, select = -label)
set.seed(42)
fit <- K_midranks(x, k = 3, dist = 1)
fit$size
fit$center
fit$converged
table(cluster = fit$cluster, truth = cauchy_k3$label)
Acc_label(fit$cluster, cauchy_k3$label)
```

Centers are in rank units. The default Manhattan variant uses coordinate-wise
medians. Set `dist = 2` for squared Euclidean dissimilarities and mean
centers, or `dist = Inf` for maximum absolute differences and coordinate-wise
maximum centers. Tied measurements receive the same rank, equal to one plus the
number of strictly larger observations in that feature. For example,
`c(9, 7, 7, 2)` receives ranks `c(1, 2, 2, 4)`. Constant features contribute
zero distance.

## Explore the number of clusters

```r
set.seed(42)
e <- elbow_rk(x, kmax = 6, verbose = FALSE, plot = TRUE)
e$objective
e$ebi
c(first = e$best.K, second = e$second.K)
```

This example fits counts 1 through 7 and considers elbows at 2 through 6.
The dashed line marks the selected elbow. The one-cluster fit supplies the baseline for scoring a possible elbow at 2.
A large ratio indicates that adding a cluster yielded a larger improvement
than adding the following one. Inspect the objective curve and convergence
flags alongside the ranking. See `?elbow_rk` for the treatment of flat or
increasing objectives.

## Reference

Larry Goldstein and Nathakhun Wiroonsri (2026). *K-Midranks: A robust clustering algorithm for heavy-tailed data with applications to spectroscopy and finance*. arXiv:2611.30953.
