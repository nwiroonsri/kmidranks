library(Kmidranks)

x <- matrix(1:4, ncol = 1)
for (distance in c(1, 2, Inf)) {
  fit <- K_midranks(x, k = 1, dist = distance)
  expected <- if (distance == 1) 4 else if (distance == 2) 5 else 6
  stopifnot(all(fit$cluster == 1L), fit$converged, fit$iter == 0L,
            fit$tot.within == expected)
}

set.seed(42)
e <- elbow_rk(x, kmax = 2, dist = 2, verbose = FALSE)
w <- e$objective[, "within"]
stopifnot(all(e$objective[, "K"] == 1:3), w[1] == 5, e$best.K == 2,
          is.na(e$second.K),
          e$ebi[1, "ratio"] == (w[1] - w[2]) / (w[2] - w[3]))
