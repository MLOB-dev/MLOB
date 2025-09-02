#' Jackknife Resampling Functions for MultiLevelOptimalBayes
#'
#' This file contains the jackknife resampling functions
#' for the MultiLevelOptimalBayes package.
#'
#' @keywords internal

#' Estimate Bayesian parameters with individual-level jackknife standard errors
#'
#' @param data List containing data for jackknife estimation
#' @return List with jackknife estimation results
#' @keywords internal
estimate_Bay_CV_SE_jackknife_individual <- function(data) {
 
  # Original estimator calculation (to get initial results)
  original_result <- estimate_Bay_CV(data)

  # Initialize parameters
  n <- data$n
  J <- data$k
  r <- min(n,50) # set up the number of replications for jackknife, r from 1 to n^J

  # Initialize vectors to store jackknife estimates for individual deletion
  jackknife_beta_b_Bay_individual <- numeric(r)
  jackknife_beta_b_Bay_ML_individual <- numeric(r)
  jackknife_beta_b_ML_individual <- numeric(r)
  
  # only if gamma-covariates present
  has_gamma <- ("kc" %in% names(data) && data$kc > 0)
  if (has_gamma) {
    p_gamma <- data$kc
    jackknife_gamma_individual <- matrix(NA, nrow = r, ncol = p_gamma)
  }

  # Jackknife by deleting one individual from each group simultaneously
  for (i in seq_len(r)) {
    # Generate J random indices (one for each group) from 1 to n
    random_indices <- sample(1:n, J, replace = TRUE)

    # Initialize a list to store indices to keep for each group
    indices_to_keep <- vector("list", J)

    for (j in 1:J) {
      indices_to_keep[[j]] <- (n * (j - 1) + 1):(n * j)
      indices_to_keep[[j]] <- indices_to_keep[[j]][-random_indices[j]]
    }

    # Create new indices after deleting one individual from each group
    indices <- unlist(indices_to_keep)

    # Create a new dataset excluding the selected individuals
    data_jackknife <- data
    data_jackknife$x <- data$x[indices]
    data_jackknife$y <- data$y[indices]
    
    if ("kc" %in% names(data)){
      if (data$kc>1){
        data_jackknife$C <- data$C[indices, , drop = FALSE]
      } else {
        data_jackknife$C <- data$C[indices, drop = FALSE]
      }
    }
    data_jackknife$n <- n - 1
    data_jackknife$kn <- (n-1)*J

    # Recalculate estimators without the selected individuals
    result_jackknife <- estimate_Bay_CV(data_jackknife)

    # Store jackknife estimates
    jackknife_beta_b_Bay_individual[i] <- result_jackknife$beta_b_Bay
    jackknife_beta_b_Bay_ML_individual[i] <- result_jackknife$beta_b_Bay_ML
    jackknife_beta_b_ML_individual[i] <- result_jackknife$beta_b_ML
    
    if (has_gamma) {
      # assume res_j$gamma is a length-kc numeric vector
      jackknife_gamma_individual[i, ] <- result_jackknife$gamma
    }
    
  }

  # Calculate jackknife means for individual deletion
  jackknife_mean_beta_b_Bay_individual <- mean(jackknife_beta_b_Bay_individual)
  jackknife_mean_beta_b_Bay_ML_individual <- mean(jackknife_beta_b_Bay_ML_individual)
  jackknife_mean_beta_b_ML_individual <- mean(jackknife_beta_b_ML_individual)
  if (has_gamma) {
    jackknife_mean_gamma_individual <- colMeans(jackknife_gamma_individual)
  }

  # Calculate jackknife standard errors for individual deletion
  
  SE <- function(vals, m) {
    sqrt((n - 1) / n * sum((vals - m)^2))
  }
  
  SE_beta_Bay_jackknife_individual       <- SE(jackknife_beta_b_Bay_individual, jackknife_mean_beta_b_Bay_individual)
  SE_beta_Bay_ML_jackknife_individual    <- SE(jackknife_beta_b_Bay_ML_individual, jackknife_mean_beta_b_Bay_ML_individual)
  SE_beta_ML_jackknife_individual        <- SE(jackknife_beta_b_ML_individual, jackknife_mean_beta_b_ML_individual)
  
  if (has_gamma) {
    
    SE_gamma_jackknife_individual <- sapply(
      seq_len(p_gamma),
      function(j) {
        SE(jackknife_gamma_individual[, j], m = jackknife_mean_gamma_individual[j])
      }
    )
    
    names(SE_gamma_jackknife_individual) <- paste0("SE_gamma_", seq_len(p_gamma))
  }
  
  # Output updated results with jackknife standard errors
  
  # 7. Return everything
  out <- list(
    beta_b_Bay                          = original_result$beta_b_Bay,
    beta_b_Bay_ML                       = original_result$beta_b_Bay_ML,
    beta_b_ML                           = original_result$beta_b_ML,
    SE_beta_Bay_jackknife_individual    = SE_beta_Bay_jackknife_individual,
    SE_beta_Bay_ML_jackknife_individual = SE_beta_Bay_ML_jackknife_individual,
    SE_beta_ML_jackknife_individual     = SE_beta_ML_jackknife_individual
  )
  
  if (has_gamma) {
    out$gamma                           <- original_result$gamma
    out$SE_gamma_jackknife_individual   <- SE_gamma_jackknife_individual
  }
  
  out
  
}
