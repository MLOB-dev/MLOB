# Bayesian Estimation Functions for MultiLevelOptimalBayes
#
# This file contains the Bayesian estimation functions
# for the MultiLevelOptimalBayes package.

#' Moments of a ratio of independent gamma-approximated variables
#'
#' Moments of \eqn{Y / X} for independent \eqn{Y} (mean `E_y`, variance `V_y`)
#' and \eqn{X \sim Gamma(k_B, \theta_B)} with mean `D = k_B theta_B` and
#' variance `v = k_B theta_B^2` (Dashuk et al. 2025a, Eqs. 70-78 and 98-101):
#' \eqn{E[1/X] = D / (D^2 - v)} and \eqn{E[1/X^2] = D^2 / ((D^2 - v)(D^2 - 2v))}.
#' Both exist only for \eqn{k_B > 2}, i.e. `D > 0` and `D^2 > 2 v`; `NA` otherwise.
#'
#' @param E_y,V_y mean and variance of the numerator.
#' @param D,v mean and variance of the denominator (vectorised).
#' @return List with the first (`m1`) and second (`m2`) raw moments.
#' @keywords internal
ratio_moments <- function(E_y, V_y, D, v) {
  ok <- D > 0 & D^2 > 2 * v
  m1 <- ifelse(ok, E_y * D / (D^2 - v), NA_real_)
  m2 <- ifelse(ok, (V_y + E_y^2) * D^2 / ((D^2 - v) * (D^2 - 2 * v)), NA_real_)
  list(m1 = m1, m2 = m2)
}

#' Standard error of a ratio of independent gamma-approximated variables
#'
#' Exact standard deviation of \eqn{Y / X} under the gamma approximation
#' (Dashuk et al. 2025a, Eqs. 102-103). Where it does not exist
#' (\eqn{k_B \le 2}, heavy tails), the first-order (delta method)
#' approximation is returned instead.
#'
#' @inheritParams ratio_moments
#' @return Standard error.
#' @keywords internal
ratio_se <- function(E_y, V_y, D, v) {
  mom <- ratio_moments(E_y, V_y, D, v)
  if (is.na(mom$m1)) {
    return(sqrt(V_y / D^2 + E_y^2 * v / D^4))
  }
  sqrt(max(mom$m2 - mom$m1^2, 0))
}

#' Estimate Bayesian parameters with control variables
#'
#' @param data List containing data for Bayesian estimation
#' @return List with Bayesian estimation results
#' @keywords internal
estimate_Bay_CV <- function(data) {

  J <- data$k
  n <- data$n

  # calculate covariance matrices:

  #cov(x):
  x <- matrix(data$x, n, J)

  av_x <- mean(x)
  av_x_j <- colMeans(x)

  SSA <- n * sum(av_x_j^2) - n * J * (av_x^2)
  SSD <- sum(data$x^2) - n * sum(av_x_j^2)

  MSA <- SSA / (J - 1)
  MSD <- SSD / ((n - 1) * J)

  tau_x2 <- (MSA - MSD) / n
  sigma_x2 <- MSD

  #Multiple control variables part
  # checks the existence of control variables and their number is not 0

  #Find 3 phi's for every of kc control variables from equation
  # C = phi1 + phi2*X_b + phi3*X_w + eps. Our aim is phi2 and phi3

  if ("kc" %in% names(data) && !is.null(data$kc) && data$kc != 0) {
    phi <- matrix(0, 3, data$kc)

    C <- list()

    for (i in 1:data$kc) {
        C[[i]] <- list()

        # 2 cases for number of control variables
        if (data$kc>1){
          # Reshape data$C[,i] into a matrix with n rows and J columns
          C[[i]]$C <- matrix(data$C[, i], nrow = n, ncol = J)
        } else{
          C[[i]]$C <- matrix(data$C, nrow = n, ncol = J)
        }

        # Compute the average of the entire C matrix
        C[[i]]$av_C <- mean(C[[i]]$C)

        # Compute the average of each column in the C matrix
        C[[i]]$av_C_j <- colMeans(C[[i]]$C)

        # Compute SPA_C and SPD_C
        C[[i]]$SPA_C <- n * sum(C[[i]]$av_C_j * av_x_j) - n * J * (av_x * C[[i]]$av_C)
        # 2 cases for number of control variables
        if (data$kc>1){
          C[[i]]$SPD_C <- sum(data$x * data$C[, i]) - n * sum(C[[i]]$av_C_j * av_x_j)
        }else{
          C[[i]]$SPD_C <- sum(data$x * data$C) - n * sum(C[[i]]$av_C_j * av_x_j)
        }
        # Compute MPA_C and MPD_C
        C[[i]]$MPA_C <- C[[i]]$SPA_C / (J - 1)
        C[[i]]$MPD_C <- C[[i]]$SPD_C / ((n - 1) * J)

        # Compute tau_C_x and sigma_C_x
        C[[i]]$tau_C_x <- (C[[i]]$MPA_C - C[[i]]$MPD_C) / n
        C[[i]]$sigma_C_x <- C[[i]]$MPD_C

        # Initialize phi as a vector of zeros
        C[[i]]$phi <- numeric(3)

        # Compute phi values
        C[[i]]$phi[2] <- C[[i]]$tau_C_x / tau_x2
        C[[i]]$phi[3] <- C[[i]]$sigma_C_x / sigma_x2

        # Store the phi vector in the phi matrix
        phi[, i] <- C[[i]]$phi
      }

    C <- data$C
    x_b <- matrix(rep(av_x_j, each = n), nrow = n * J, ncol = 1)

    x_w <- as.vector(matrix(data$x, nrow = n * J, ncol = 1)) - x_b
    C_tilde <- C - cbind(x_b, x_w) %*% phi[2:nrow(phi), ]
    C_tilde_mat <- cbind(rep(1, data$kn), C_tilde)

    #Estimate gamma from equation y = gamma1 + gamma*C_tilde + eps. Our
    #interest is vestor of parameters gamma without intersection gamma1
    gamma_all <- solve(t(C_tilde_mat) %*% C_tilde_mat) %*% t(C_tilde_mat) %*% data$y

    # OLS standard errors of the same regression (intercept and all covariates)
    res_gamma <- data$y - C_tilde_mat %*% gamma_all
    var_res_gamma <- sum(res_gamma^2) / (data$kn - data$kc - 1)
    SE_gamma <- sqrt(diag(var_res_gamma * solve(crossprod(C_tilde_mat))))[-1]

    gamma <- gamma_all[-1]
    # Find y_tilde - the difference between real and estimated
    # y (without intersection)
    if (data$kc>1){
      y_tilde <- data$y - C %*% gamma
    }else{
      y_tilde <- data$y - C * gamma
    }
    data$y <- y_tilde
  }

  #cov(x,y)
  y <- matrix(data$y, n, J)

  av_y <- mean(y)
  av_y_j <- colMeans(y)

  SPA <- n * sum(av_y_j * av_x_j) - n * J * (av_x * av_y)
  SPD <- sum(data$x * data$y) - n * sum(av_y_j * av_x_j)

  MPA <- SPA / (J - 1)
  MPD <- SPD / ((n - 1) * J)

  tau_yx <- (MPA - MPD) / n
  sigma_yx <- MPD

  #cov(y):
  SOA <- n * sum(av_y_j^2) - n * J * (av_y^2)
  SOD <- sum(data$y^2) - n * sum(av_y_j^2)

  MOA <- SOA / (J - 1)
  MOD <- SOD / ((n - 1) * J)

  tau_y2 <- (MOA - MOD) / n
  sigma_y2 <- MOD

  # Mean and variance of tau_x2 = (MSA - MSD) / n and tau_yx = (MPA - MPD) / n
  # under the two-level normal model, with the estimated components plugged in.
  # The between- and within-group mean-square matrices are independent Wishart
  # with J - 1 and J(n - 1) degrees of freedom; their expectations are estimated
  # by the observed mean squares. These are the exact moments of the quadratic
  # forms that Dashuk et al. (2025a, Eqs. 104-105) approximate by gamma
  # distributions (Welch-Satterthwaite: k = E^2 / V, theta = V / E).
  df_A <- J - 1
  df_D <- J * (n - 1)

  E_x  <- tau_x2
  V_x  <- 2 * (MSA^2 / df_A + MSD^2 / df_D) / n^2
  E_yx <- tau_yx
  V_yx <- ((MPA^2 + MSA * MOA) / df_A + (MPD^2 + MSD * MOD) / df_D) / n^2

  beta_b_ML <- tau_yx / tau_x2

  # MSE of the ML estimator, evaluated at the true beta_b when supplied (simulations)
  beta_b <- if (is.null(data$b2)) beta_b_ML else data$b2
  mom_ML <- ratio_moments(E_yx, V_yx, E_x, V_x)
  MSE_add <- mom_ML$m2 - 2 * beta_b * mom_ML$m1 + beta_b^2

  # optimize MSE with beta_b_ML within a 5 SD search region from the distribution of tau_x2

  d_search <- 5
  radius   <- sqrt(V_x)

  # Grid bounds are relative to the estimated total variance of x, so that the
  # estimator does not depend on the units of x. For x with unit total variance
  # the floor equals 0.01 (a between-group share of 1%).
  s2_x        <- tau_x2 + sigma_x2
  Tau02_floor <- 0.01 * s2_x

  # Attempt the usual endpoints
  Tau02_min <- max(Tau02_floor, tau_x2 - d_search * radius)
  Tau02_max <- max(Tau02_floor, tau_x2 + d_search * radius)

  # Fallback to a fixed window (between-group share of 1% to 100%) if anything's off
  if (!all(is.finite(c(Tau02_min, Tau02_max))) || Tau02_min >= Tau02_max) {
    Tau02_min <- Tau02_floor
    Tau02_max <- s2_x
    warning('Please check grouping, as it does not explain the variation and may be meaningless.')
  }

  # Build the sequence
  Tau02_seq <- seq(
    from       = Tau02_min,
    to         = Tau02_max,
    length.out = 100
  )

  W_seq <- seq(0.01, 1, by = 0.01)

  grid <- expand.grid(Tau02 = Tau02_seq, W = W_seq)

  Tau02 <- grid$Tau02

  W <- grid$W

  # Denominator X = (1 - w) tau02 + w tau_x2: mean D, variance v
  D <- (1 - W) * Tau02 + W * E_x
  v <- W^2 * V_x
  mom <- ratio_moments(E_yx, V_yx, D, v)

  # MSE(w, tau02) = E[beta^2] - 2 b E[beta] + b^2 (Dashuk et al. 2025a, Eqs. 78 and 101);
  # grid points where the moments do not exist are excluded
  mse_grid <- function(b) {
    mse <- mom$m2 - 2 * b * mom$m1 + b^2
    mse[is.na(mse)] <- Inf
    mse
  }

  MSE_ML <- mse_grid(beta_b_ML)  # with the ML estimate plugged in for beta_b
  MSE    <- mse_grid(data$b_b)   # with the value of beta_b supplied in data$b_b

  if (!any(is.finite(MSE)) || !any(is.finite(MSE_ML))) {
    stop("The MSE of the Bayesian estimator cannot be evaluated on the search grid.")
  }

  # use MSE to find beta_b_Bay with smallest MSE:

  ind_MSEbetabBay <- which.min(MSE)
  MSE_beta_b_Bay <- MSE[ind_MSEbetabBay]

  w_opt <- W[ind_MSEbetabBay]
  tau02_opt <- Tau02[ind_MSEbetabBay]

  beta_b_Bay <- tau_yx/((1-w_opt)*tau02_opt + w_opt*tau_x2)

  # use MSE_ML to find beta_b_Bay_ML with smallest MSE:

  ind_MSEbetabBayML <- which.min(MSE_ML)
  MSE_beta_b_Bay_ML <- MSE_ML[ind_MSEbetabBayML]

  w_opt_ML <- W[ind_MSEbetabBayML]
  tau02_opt_ML <- Tau02[ind_MSEbetabBayML]

  beta_b_Bay_ML <- tau_yx/((1-w_opt_ML)*tau02_opt_ML + w_opt_ML*tau_x2)

  # Standard errors of beta_b_ML, beta_b_Bay and beta_b_Bay_ML from their distributions

  SE_beta_ML <- ratio_se(E_yx, V_yx, E_x, V_x)

  SE_beta_Bay <- ratio_se(E_yx, V_yx, (1 - w_opt) * tau02_opt + w_opt * E_x,
                          w_opt^2 * V_x)

  SE_beta_Bay_ML <- ratio_se(E_yx, V_yx, (1 - w_opt_ML) * tau02_opt_ML + w_opt_ML * E_x,
                             w_opt_ML^2 * V_x)

  Bay <- list()
  Bay$tau_x2 <- tau_x2
  Bay$sigma_x2 <- sigma_x2
  Bay$tau_yx <- tau_yx
  Bay$sigma_yx <- sigma_yx
  Bay$tau_y2 <- tau_y2
  Bay$sigma_y2 <- sigma_y2
  Bay$E_x <- E_x
  Bay$V_x <- V_x
  Bay$E_yx <- E_yx
  Bay$V_yx <- V_yx
  Bay$MSE <- MSE
  Bay$MSE_ML <- MSE_ML
  Bay$W <- W
  Bay$Tau02 <- Tau02
  Bay$w_opt <- w_opt
  Bay$tau02_opt <- tau02_opt
  Bay$tau02_at_edge <- tau02_opt %in% c(Tau02_min, Tau02_max)
  Bay$w_opt_ML <- w_opt_ML
  Bay$tau02_opt_ML <- tau02_opt_ML
  Bay$beta_b_ML <- beta_b_ML
  Bay$beta_b_Bay <- beta_b_Bay
  Bay$beta_b_Bay_ML <- beta_b_Bay_ML

  if (exists("gamma") && !is.function(gamma)){
    Bay$gamma <- gamma
  }

  Bay$MSE_beta_b_Bay <- MSE_beta_b_Bay
  Bay$MSE_beta_b_Bay_ML <- MSE_beta_b_Bay_ML
  Bay$MSE_add <- MSE_add
  Bay$SE_beta_ML = SE_beta_ML
  Bay$SE_beta_Bay = SE_beta_Bay
  Bay$SE_beta_Bay_ML = SE_beta_Bay_ML

  if (exists("SE_gamma")){
    Bay$SE_gamma = SE_gamma
  }

  return(Bay)
}
