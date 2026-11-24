
K_midranks <- function(x, k, dist = 1, nstart = 20, mr = 1000) {
  
  # --------------------------------------------------
  # 1. Input validation
  # --------------------------------------------------
  
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a numeric matrix or data frame.")
  }
  
  if (!all(vapply(as.data.frame(x), is.numeric, logical(1)))) {
    stop("All columns of 'x' must be numeric.")
  }
  
  x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  
  if (n < 2 || p < 1) {
    stop("'x' must have at least 2 rows and 1 column.")
  }
  
  if (anyNA(x) || any(!is.finite(x))) {
    stop("'x' must contain only finite values.")
  }
  
  check_integer <- function(value, name, lower) {
    if (!is.numeric(value) || length(value) != 1 ||
        is.na(value) || !is.finite(value) ||
        value != floor(value) || value < lower) {
      stop(sprintf("'%s' must be an integer >= %d.",
                   name, lower))
    }
  }
  
  check_integer(k, "k", 1)
  check_integer(nstart, "nstart", 1)
  check_integer(mr, "mr", 1)
  
  if (k > n) {
    stop("'k' cannot exceed the number of observations.")
  }
  
  if (!is.numeric(dist) || length(dist) != 1 ||
      is.na(dist) || !(dist %in% c(1, 2, Inf))) {
    stop("'dist' must be 1, 2, or Inf.")
  }
  
  if (any(vapply(seq_len(p),
                 function(j) length(unique(x[, j])) == 1,
                 logical(1)))) {
    warning("Some variables are constant.")
  }
  
  # --------------------------------------------------
  # 2. Rank transformation
  # --------------------------------------------------
  
  rk <- vapply(
    seq_len(p),
    function(j) rank(-x[, j], ties.method = "min"),
    numeric(n)
  )
  
  # --------------------------------------------------
  # 3. Helper functions
  # --------------------------------------------------
  
  update_centers <- function(gg) {
    
    cc <- matrix(0, nrow = k, ncol = p)
    
    for (i in seq_len(k)) {
      
      cluster_data <- rk[gg == i, , drop = FALSE]
      
      if (dist == 1) {
        cc[i, ] <- apply(cluster_data, 2, median)
        
      } else if (dist == 2) {
        cc[i, ] <- colMeans(cluster_data)
        
      } else {
        # Preserve original L-infinity center update
        cc[i, ] <- apply(cluster_data, 2, max)
      }
    }
    
    cc
  }
  
  calculate_distances <- function(cc) {
    
    aa <- matrix(0, nrow = n, ncol = k)
    
    for (i in seq_len(k)) {
      
      dif <- abs(sweep(rk, 2, cc[i, ], "-"))
      
      if (dist == 1) {
        aa[, i] <- rowSums(dif)
        
      } else if (dist == 2) {
        # Squared Euclidean distance
        aa[, i] <- rowSums(dif^2)
        
      } else {
        # Maximum absolute coordinate difference
        aa[, i] <- apply(dif, 1, max)
      }
    }
    
    aa
  }
  
  # Preserve original empty-cluster repair strategy
  repair_empty_clusters <- function(gg) {
    
    counts <- tabulate(gg, nbins = k)
    empty_clusters <- which(counts == 0)
    
    if (length(empty_clusters) == 0) {
      return(gg)
    }
    
    sss <- setdiff(seq_len(k), empty_clusters)
    
    for (empty in empty_clusters) {
      
      eligible <- sss[counts[sss] >= 2]
      
      if (length(eligible) == 0) {
        stop("Unable to repair empty clusters.")
      }
      
      donor <- eligible[1]
      
      # Move the first observation from the donor
      idx <- which(gg == donor)[1]
      gg[idx] <- empty
      
      counts[donor] <- counts[donor] - 1
      counts[empty] <- 1
      
      # Rotate donor to the end
      sss <- c(sss[sss != donor], donor)
    }
    
    gg
  }
  
  # --------------------------------------------------
  # 4. Multiple random starts
  # --------------------------------------------------
  
  if (k == 1L) {
    cluster <- rep.int(1L, n)
    center <- update_centers(cluster)
    within <- as.vector(calculate_distances(center))
    return(list(cluster = cluster, center = center, size = n, within = within,
                tot.within = sum(within), iter = 0L, converged = TRUE))
  }

  best_objective <- Inf
  best_result <- NULL
  
  for (ns in seq_len(nstart)) {
    
    # Random initialization with no empty clusters
    gg <- sample(rep(seq_len(k), length.out = n))
    
    cc <- update_centers(gg)
    converged <- FALSE
    
    for (r in seq_len(mr)) {
      
      # Assignment step
      aa <- calculate_distances(cc)
      gg_new <- max.col(-aa, ties.method = "first")
      
      # Repair empty clusters
      gg_new <- repair_empty_clusters(gg_new)
      
      # Update cluster centers
      cc_new <- update_centers(gg_new)
      
      # Calculate within-cluster distances
      aa_new <- calculate_distances(cc_new)
      wcd <- aa_new[cbind(seq_len(n), gg_new)]
      
      # Check convergence
      if (all(gg_new == gg)) {
        converged <- TRUE
        gg <- gg_new
        cc <- cc_new
        break
      }
      
      gg <- gg_new
      cc <- cc_new
    }
    
    # ------------------------------------------------
    # 5. Save the best solution
    # ------------------------------------------------
    
    objective <- sum(wcd)
    
    if (objective < best_objective) {
      
      best_objective <- objective
      
      best_result <- list(
        cluster = gg,
        center = cc,
        size = tabulate(gg, nbins = k),
        within = wcd,
        tot.within = objective,
        iter = r,
        converged = converged
      )
    }
  }
  
  if (!best_result$converged) {
    warning(
      "The best solution did not converge within 'mr' iterations."
    )
  }
  
  return(best_result)
}

