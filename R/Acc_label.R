
Acc_label <- function(kmc, y) {
  
  # --------------------------------------------------
  # 1. Input validation
  # --------------------------------------------------
  
  if (!is.atomic(kmc) || !is.atomic(y) ||
      !is.null(dim(kmc)) || !is.null(dim(y))) {
    stop("'kmc' and 'y' must be vectors.")
  }
  
  if (length(kmc) != length(y)) {
    stop("'kmc' and 'y' must have the same length.")
  }
  
  if (length(kmc) == 0) {
    stop("'kmc' and 'y' must not be empty.")
  }
  
  if (anyNA(kmc) || anyNA(y)) {
    stop("'kmc' and 'y' must not contain missing values.")
  }
  
  # --------------------------------------------------
  # 2. Construct contingency table
  # --------------------------------------------------
  
  tab <- table(factor(kmc), factor(y))
  
  # --------------------------------------------------
  # 3. Optimal label matching
  # --------------------------------------------------
  
  # solve_LSAP requires rows <= columns.
  # Pad the table with zeros if necessary.
  
  nr <- nrow(tab)
  nc <- ncol(tab)
  
  if (nr > nc) {
    tab <- cbind(tab, matrix(0, nr, nr - nc))
  }
  
  perm <- clue::solve_LSAP(tab, maximum = TRUE)
  
  # --------------------------------------------------
  # 4. Calculate accuracy
  # --------------------------------------------------
  
  matched <- sum(tab[cbind(seq_len(nr), as.integer(perm))])
  
  accuracy <- matched / length(y)
  
  return(accuracy)
}