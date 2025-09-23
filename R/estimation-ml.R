#' ML Estimation Functions for MultiLevelOptimalBayes
#'
#' This file contains the Maximum Likelihood estimation functions
#' for the MultiLevelOptimalBayes package.
#'
#' @keywords internal

#' Estimate ML parameters with control variables
#'
#' @param data List containing data for ML estimation
#' @return List with ML estimation results
#' @keywords internal
estimate_ML_CV <- function(data) {

  data$j <- data$k
  # Reshape the data.x into an n x j matrix
  x <- matrix(data$x, nrow = data$n, ncol = data$k)

  av_x <- mean(x)
  av_x_j <- colMeans(x)

  SSA <- data$n * sum(av_x_j^2) - data$n * data$k * (av_x^2)
  SSD <- sum(data$x * data$x) - data$n * sum(av_x_j * av_x_j)

  MSA <- SSA / (data$j - 1)
  MSD <- SSD / ((data$n - 1) * data$j)

  tau_x2 <- (MSA - MSD) / data$n
  sigma_x2 <- MSD

  # Reshape the data.y into an n x j matrix
  y <- matrix(data$y, nrow = data$n, ncol = data$j)

  av_y <- mean(y)
  av_y_j <- colMeans(y)

  SPA <- data$n * sum(av_y_j * av_x_j) - data$n * data$j * (av_x * av_y)
  SPD <- sum(data$x *data$y) - data$n * sum(av_y_j * av_x_j)

  MPA <- SPA / (data$j - 1)
  MPD <- SPD / ((data$n - 1) * data$j)

  tau_yx <- (MPA - MPD) / data$n
  sigma_yx <- MPD

  # ML estimation of beta_b
  beta_b_ML <- tau_yx / tau_x2

  # Multiple control variables part

  if (!is.null(data$kc) && data$kc != 0) {

    beta_b_ML_tilde <- NULL
    data$C = matrix(data$C)

    phi <- matrix(0, nrow = 3, ncol = data$kc)
    data$C = matrix(data$C, nrow = data$k*data$n, ncol = data$kc)
    for (i in 1:data$kc) {

      if (data$kc>1){
        # Reshape data$C[,i] into a matrix with n rows and J columns
        C_i_C <- matrix(data$C[, i], nrow = data$n, ncol = data$k)
      }else{
        C_i_C <- matrix(data$C, nrow = data$n, ncol = data$k)
      }

      av_C <- mean(C_i_C)
      av_C_j <- colMeans(C_i_C)

      SPA_C <- data$n * sum(av_C_j * av_x_j) - data$n * data$j * (av_x * av_C)
      SPD_C <- sum(x * data$C[, i]) - data$n * sum(av_C_j * av_x_j)

      MPA_C <- SPA_C / (data$j - 1)
      MPD_C <- SPD_C / ((data$n - 1) * data$j)

      tau_C_x <- (MPA_C - MPD_C) / data$n
      sigma_C_x <- MPD_C

      phi[2, i] <- tau_C_x / tau_x2
      phi[3, i] <- sigma_C_x / sigma_x2
    }

    C <- data$C
    x_b <- rep(av_x_j, each = data$n)
    x_w <- as.vector(data$x) - x_b
    C_tilde <- C - cbind(x_b, x_w) %*% (phi[2:3, ])

    C_tilde_mat <- cbind(1, C_tilde)

    gamma_all <- solve(t(C_tilde_mat) %*% C_tilde_mat) %*% t(C_tilde_mat) %*% data$y
    gamma <- gamma_all[-1]

    y_tilde <- data$y - C %*% gamma
    y_tilde <- matrix(y_tilde, nrow = data$n, ncol = data$j)

    av_y_tilde <- mean(y_tilde)
    av_y_tilde_j <- colMeans(y_tilde)

    SPA_tilde <- data$n * sum(av_y_tilde_j * av_x_j) - data$n * data$j * (av_x * av_y_tilde)
    SPD_tilde <- sum(as.vector(data$x) * as.vector(y_tilde)) - data$n * sum(av_y_tilde_j * av_x_j)

    MPA_tilde <- SPA_tilde / (data$j - 1)
    MPD_tilde <- SPD_tilde / ((data$n - 1) * data$j)

    tau_y_tilde_x <- (MPA_tilde - MPD_tilde) / data$n

    # ML estimation of beta_b for the model with control variables
    beta_b_ML_tilde <- tau_y_tilde_x / tau_x2

    return(list(beta_b_ML_CV = beta_b_ML_tilde, beta_b_ML_no_CV = beta_b_ML, 
                gamma = gamma, tau_x2 = tau_x2, sigma_x2 = sigma_x2, 
                tau_yx = tau_yx, sigma_yx = sigma_yx))

  } else {

    return(list(beta_b_ML_CV = beta_b_ML, beta_b_ML_no_CV = beta_b_ML, 
                tau_x2 = tau_x2, sigma_x2 = sigma_x2, tau_yx = tau_yx, 
                sigma_yx = sigma_yx))
  }
}
