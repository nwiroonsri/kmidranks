
elbow_rk <- function(x, dist = 1, kmax = 10,
                     nstart = 5, mr = 1000,
                     verbose = TRUE, plot = FALSE) {
  
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
  
  if (nrow(x) < 3 || ncol(x) < 1) {
    stop("'x' must have at least 3 rows and 1 column.")
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
  
  check_integer(kmax, "kmax", 2)
  check_integer(nstart, "nstart", 1)
  check_integer(mr, "mr", 1)
  
  if (kmax + 1 > nrow(x)) {
    stop("'kmax + 1' cannot exceed the number of observations.")
  }
  
  if (!is.numeric(dist) || length(dist) != 1 ||
      is.na(dist) || !(dist %in% c(1, 2, Inf))) {
    stop("'dist' must be 1, 2, or Inf.")
  }
  
  if (!is.logical(verbose) || length(verbose) != 1 ||
      is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.")
  }
  
  # --------------------------------------------------
  # 2. Calculate within-cluster objectives
  # --------------------------------------------------
  
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("'plot' must be TRUE or FALSE.")
  }

  start_time <- proc.time()[["elapsed"]]
  
  ks <- seq_len(kmax + 1)
  eb <- numeric(length(ks))
  converged <- logical(length(ks))
  
  for (i in seq_along(ks)) {
    
    k <- ks[i]
    
    fit <- K_midranks(
      x,
      k = k,
      dist = dist,
      nstart = nstart,
      mr = mr
    )
    
    eb[i] <- sum(fit$within)
    converged[i] <- fit$converged
    
    if (verbose) {
      pct <- round(100 * i / length(ks))
      cat(sprintf("\rElbow calculation: %3d%% (K = %d)",
                  pct, k))
      flush.console()
    }
  }
  
  if (verbose) cat("\n")
  
  if (any(!converged)) {
    warning("Some clustering solutions did not converge.")
  }
  
  # --------------------------------------------------
  # 3. Calculate slopes
  # --------------------------------------------------
  
  slope <- eb[-length(eb)] - eb[-1]
  
  # Ratio associated with candidate K = 2,...,kmax
  candidate_k <- 2:kmax
  
  numerator <- slope[-length(slope)]
  denominator <- slope[-1]
  
  ratio <- rep(NA_real_, length(candidate_k))
  
  # Avoid division by zero and undefined ratios
  valid <- is.finite(numerator) &
    is.finite(denominator) &
    numerator >= 0 &
    denominator > 0
  
  ratio[valid] <- numerator[valid] / denominator[valid]
  
  if (any(!valid)) {
    warning(
      "Some elbow ratios are undefined or invalid; ",
      "their values are returned as NA."
    )
  }
  
  ebp <- cbind(K = candidate_k, ratio = ratio)
  
  # --------------------------------------------------
  # 4. Identify the two largest elbow ratios
  # --------------------------------------------------
  
  valid_idx <- which(is.finite(ratio))
  
  if (length(valid_idx) == 0) {
    best_k <- NA_integer_
    second_k <- NA_integer_
    
    warning("No valid elbow ratio was found.")
    
  } else {
    
    ord <- valid_idx[order(ratio[valid_idx],
                           decreasing = TRUE)]
    
    best_k <- candidate_k[ord[1]]
    
    second_k <- if (length(ord) >= 2) {
      candidate_k[ord[2]]
    } else {
      NA_integer_
    }
  }
  
  # --------------------------------------------------
  # 5. Return results
  # --------------------------------------------------
  
  elapsed <- proc.time()[["elapsed"]] - start_time
  
  my_list <- list(
    best.K = best_k,
    second.K = second_k,
    ebv = cbind(K = ks[-length(ks)], slope = slope),
    ebi = ebp,
    objective = cbind(K = ks, within = eb),
    converged = cbind(K = ks, converged = converged),
    elapsed = elapsed
  )
  
  if (plot) {
    graphics::plot(ks, eb, type = "b", xlab = "Number of clusters (K)",
                   ylab = "Total within-cluster dissimilarity")
    if (!is.na(best_k)) {
      graphics::abline(v = best_k, lty = 2, col = "steelblue4")
      graphics::points(best_k, eb[match(best_k, ks)], pch = 19,
                       col = "steelblue4")
    }
  }

  return(my_list)
}
