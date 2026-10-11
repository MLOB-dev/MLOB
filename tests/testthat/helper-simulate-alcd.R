# Test fixture: R port of the published ALCD data-generating process
# (Dashuk et al. 2024/2025b supplementary material).
#
# Ports two MATLAB functions 1:1, generalized to any number of control
# variables (kc = length(gamma)):
#   - assembleLatentCovariateData.m  -> kc = 0 branch (no control variables)
#   - ALCD_CV.m                      -> kc >= 1 branch (kc control variables)
#
# ALCD_CV.m's own comment on "kc" says: "if change kc, also change gamma, m_C,
# ICC_C (all vectors), cov_mat, cov_mat_b (matrices)". This port follows that
# literally: kc is *inferred* from length(gamma) rather than passed
# separately, and gamma/m_C/ICC_C/cov_mat/cov_mat_b are ordinary arguments, so
# it works for kc = 0, 1, 3 (the two MATLAB configurations seen so far), or
# any other kc, as long as the caller supplies matrices of matching dimension.
#
# By construction (no sampling error, these are the *population* values used
# to build the data):
#   sigma_x2 = 1 - ICC_x,  tau_x2 = ICC_x
#   tau_yx   = beta_b * tau_x2,   sigma_yx = beta_w * sigma_x2
# so beta_b_ML = tau_yx / tau_x2 = beta_b exactly, in expectation, regardless
# of kc.
#
# testthat automatically sources files named "helper-*.R" before running
# tests, and does NOT install them as part of the package.

#' Draw from a multivariate normal via a Cholesky factor
#'
#' Same construction as MATLAB's `mvnrnd` (`randn(n, p) * chol(Sigma) + mu`).
#' Used instead of `MASS::mvrnorm()`, whose eigendecomposition-based draws can
#' differ across BLAS/LAPACK builds (eigenvector signs are not unique), which
#' would make seeded test data platform-dependent. Falls back to a pivoted
#' Cholesky factor for positive semidefinite (singular) `Sigma`.
rmvnorm_chol <- function(n, mu, Sigma) {
  p <- ncol(Sigma)
  R <- tryCatch(chol(Sigma), error = function(e) {
    Rp <- suppressWarnings(chol(Sigma, pivot = TRUE))
    rank <- attr(Rp, "rank")
    if (rank < p) Rp[(rank + 1):p, (rank + 1):p] <- 0
    Rp[, order(attr(Rp, "pivot")), drop = FALSE]
  })
  matrix(stats::rnorm(n * p), nrow = n) %*% R + matrix(mu, nrow = n, ncol = p, byrow = TRUE)
}

#' Simulate two-level data from the published ALCD / ALCD_CV DGP
#'
#' @param J number of groups (`k`/`j` in the MATLAB source and Zitzmann et al. 2020).
#' @param n number of observations per group.
#' @param ICC_x target intraclass correlation of the predictor `x`.
#' @param beta_b true between-group coefficient (`b2`/`beta_b` in the MATLAB source).
#' @param beta_w true within-group coefficient (`b1`/`beta_w` in the MATLAB source).
#' @param ICC_y target intraclass correlation of the response `y`. Fixed at 0.2
#'   in both MATLAB sources; exposed here as an argument for flexibility.
#' @param gamma numeric vector of true control-variable coefficients. `NULL`
#'   (default) reproduces `assembleLatentCovariateData.m` (no control
#'   variables). A vector of length `kc` reproduces `ALCD_CV.m` with `kc`
#'   control variables -- `kc` is inferred as `length(gamma)`, not passed
#'   separately, mirroring the MATLAB source's own comment that `gamma`,
#'   `m_C`, `ICC_C`, `cov_mat`, and `cov_mat_b` must all change together with `kc`.
#' @param m_C numeric vector of length `kc`, population means of the control
#'   variables. Defaults to all zeros if `gamma` is supplied but `m_C` is not.
#' @param ICC_C numeric vector of length `kc`, target ICCs of the control variables.
#' @param cov_mat a `(kc+1) x (kc+1)` total covariance matrix of `(x, C_1, ..., C_kc)`,
#'   in that variable order (matches MATLAB's `cov_mat`).
#' @param cov_mat_b the corresponding `(kc+1) x (kc+1)` between-group covariance
#'   matrix (matches MATLAB's `cov_mat_b`); its diagonal must equal
#'   `c(ICC_x, ICC_C)`. `cov_mat - cov_mat_b` is used as the within-group
#'   covariance matrix (MATLAB's `cov_mat_w`).
#'
#' @return A data frame with columns `y`, `x` (and `C1, ..., Ckc` when `kc >= 1`),
#'   and `group`, ready to be passed to `mlob(y ~ x + ..., data = ..., group = "group")`.
#'   The exact population values used to construct the data (e.g. `tau_x2`,
#'   `tau_yx`) are attached as the `"true"` attribute -- see [attr()].
simulate_alcd_cv <- function(J, n, ICC_x, beta_b, beta_w, ICC_y = 0.2,
                              gamma = NULL, m_C = NULL, ICC_C = NULL,
                              cov_mat = NULL, cov_mat_b = NULL) {

  k  <- J
  b0 <- 0
  b1 <- beta_w
  b2 <- beta_b
  kc <- length(gamma)

  # -- input validation, ported from both MATLAB sources --
  if (ICC_x < max(1 - (1 - ICC_y) / b1^2, 0)) {
    stop(sprintf("Too small value of ICC_x. The value should be larger than %f",
                  max(1 - (1 - ICC_y) / b1^2, 0)))
  }
  if (ICC_x > min(ICC_y / b2^2, 1)) {
    stop(sprintf("Too large value of ICC_x. The value should be smaller than %f",
                  min(ICC_y / b2^2, 1)))
  }
  if (k <= 0 || n <= 0) {
    stop("'J' or 'n' are negative. Consider natural numbers.")
  }
  if (k != round(k) || n != round(n)) {
    stop("'J' or 'n' are non-integers. Consider using natural numbers.")
  }

  kn  <- k * n
  g   <- rep(seq_len(k), each = n)   # group membership, kron((1:k)', ones(n,1)) in MATLAB
  m_x <- 0

  if (kc == 0) {

    # ---- ported from assembleLatentCovariateData.m ----
    var_x1 <- 1 - ICC_x                            # within variance sigma_x2
    var_x2 <- ICC_x                                # between variance tau_x2
    var_e2 <- ICC_y - b2^2 * ICC_x                 # between residual variance
    var_e1 <- (1 - ICC_y) - b1^2 * (1 - ICC_x)     # within residual variance

    if (var_e1 <= 0 || var_e2 <= 0) {
      stop("Values of ICC_x and ICC_y do not coincide (negative residual variance).")
    }

    x2 <- stats::rnorm(k,  mean = m_x, sd = sqrt(var_x2))  # between part of x, k elements
    e2 <- stats::rnorm(k,  mean = 0,   sd = sqrt(var_e2))  # between part residual, k elements
    x  <- stats::rnorm(kn, mean = 0,   sd = sqrt(var_x1)) + x2[g]  # overall x, kn elements

    y <- stats::rnorm(kn, mean = 0, sd = sqrt(var_e1)) +
      (b0 + b1 * (x - x2[g]) + b2 * x2[g] + e2[g])

    df <- data.frame(y = y, x = x, group = g)

    true <- list(J = k, n = n, kn = kn, ICC_x = ICC_x, ICC_y = ICC_y,
                 beta_b = b2, beta_w = b1,
                 tau_x2 = var_x2, sigma_x2 = var_x1,
                 tau_yx = b2 * var_x2, sigma_yx = b1 * var_x1)

  } else {

    # ---- ported from ALCD_CV.m, generalized to any kc = length(gamma) ----
    gamma <- as.numeric(gamma)
    if (is.null(m_C)) m_C <- rep(0, kc)
    m_C   <- as.numeric(m_C)
    ICC_C <- as.numeric(ICC_C)

    if (length(m_C) != kc || length(ICC_C) != kc) {
      stop("'gamma', 'm_C', and 'ICC_C' must all have the same length (= kc).")
    }
    if (is.null(cov_mat) || is.null(cov_mat_b)) {
      stop("For kc >= 1 control variables, 'cov_mat' and 'cov_mat_b' ((kc+1) x (kc+1), ",
           "ordered [x, C1, ..., Ckc]) must be supplied -- see ALCD_CV.m's own comment: ",
           "'if change kc, also change gamma, m_C, ICC_C, cov_mat, cov_mat_b'.")
    }
    if (!all(dim(cov_mat) == kc + 1) || !all(dim(cov_mat_b) == kc + 1)) {
      stop(sprintf("'cov_mat' and 'cov_mat_b' must be %d x %d (kc + 1), matching length(gamma) = %d.",
                    kc + 1, kc + 1, kc))
    }
    if (max(abs(diag(cov_mat_b) - c(ICC_x, ICC_C))) > 1e-8) {
      stop("diag(cov_mat_b) must equal c(ICC_x, ICC_C).")
    }

    cov_mat_w <- cov_mat - cov_mat_b

    if (!isSymmetric(cov_mat) || !isSymmetric(cov_mat_b) || !isSymmetric(cov_mat_w)) {
      stop("Covariance matrix is not symmetric.")
    }
    if (any(eigen(cov_mat,   only.values = TRUE)$values <= -0.00001) ||
        any(eigen(cov_mat_b, only.values = TRUE)$values <= -0.00001) ||
        any(eigen(cov_mat_w, only.values = TRUE)$values <= -0.00001)) {
      stop("Covariance matrix is not positive semidefinite.")
    }

    var_x1 <- 1 - ICC_x     # within variance sigma_x2
    var_x2 <- ICC_x         # between variance tau_x2
    var_C1 <- diag(1 - ICC_C, nrow = kc)  # within variance matrix of control variables
    var_C2 <- diag(ICC_C, nrow = kc)      # between variance matrix of control variables

    # extract covariances as vectors, column by column below the diagonal, e.g.
    # cov_xC = [cov(x,C1); cov(x,C2); cov(x,C3); cov(C1,C2); cov(C1,C3); cov(C2,C3)]
    cov_xC <- cov_mat[lower.tri(cov_mat)]      # overall covariances between x and all C's
    ICC_xC <- cov_mat_b[lower.tri(cov_mat_b)]  # between-part covariances between x and all C's

    # all pairwise products of coefficients, in the same order as cov_xC/ICC_xC
    # above (mirrors MATLAB's prod(nchoosek([b; gamma], 2)'))
    pair_products <- function(v) {
      if (length(v) < 2) return(numeric(0))
      apply(utils::combn(v, 2), 2, prod)
    }

    var_e2 <- ICC_y - b2^2 * ICC_x - sum(gamma^2 * ICC_C) -
      sum(pair_products(c(b2, gamma)) * ICC_xC)              # between residual variance
    var_e1 <- (1 - ICC_y) - b1^2 * (1 - ICC_x) - sum(gamma^2 * (1 - ICC_C)) -
      sum(pair_products(c(b1, gamma)) * (cov_xC - ICC_xC))   # within residual variance

    if (var_e1 <= 0 || var_e2 <= 0) {
      stop("Values of ICC_x, ICC_C, and ICC_y do not coincide.")
    }

    # random data generation: [x, C1, ..., Ckc] jointly, group-level then overall
    xC2 <- rmvnorm_chol(k,  mu = c(m_x, m_C),    Sigma = cov_mat_b)  # k x (kc+1)
    xC  <- rmvnorm_chol(kn, mu = rep(0, kc + 1), Sigma = cov_mat_w) + xC2[g, , drop = FALSE]

    x2 <- xC2[, 1]                       # between part of x, k elements
    x  <- xC[, 1]                        # overall x, kn elements
    C  <- xC[, -1, drop = FALSE]         # overall control variables, kn x kc
    e2 <- stats::rnorm(k, mean = 0, sd = sqrt(var_e2))  # between part residual, k elements

    y <- stats::rnorm(kn, mean = 0, sd = sqrt(var_e1)) +
      (b0 + b1 * (x - x2[g]) + b2 * x2[g] + e2[g]) + C %*% gamma

    df <- data.frame(y = as.numeric(y), x = x, C, group = g)
    colnames(df)[3:(2 + kc)] <- paste0("C", seq_len(kc))

    true <- list(J = k, n = n, kn = kn, ICC_x = ICC_x, ICC_y = ICC_y, ICC_C = ICC_C,
                 beta_b = b2, beta_w = b1, gamma = gamma,
                 tau_x2 = var_x2, sigma_x2 = var_x1,
                 tau_yx = b2 * var_x2, sigma_yx = b1 * var_x1,
                 var_C1 = var_C1, var_C2 = var_C2)
  }

  attr(df, "true") <- true
  df
}

#' Published parameter configurations for `simulate_alcd_cv()`'s covariate branch
#'
#' Returns the exact `gamma`/`m_C`/`ICC_C`/`cov_mat`/`cov_mat_b` constants
#' hardcoded in the MATLAB source `ALCD_CV.m`, for two configurations seen in
#' that file: one control variable (`kc = 1`) and three control variables
#' (`kc = 3`). Convenience only -- `simulate_alcd_cv()` itself works for any
#' `kc` given matching arguments.
#'
#' @param ICC_x the same `ICC_x` that will be passed to `simulate_alcd_cv()`
#'   (needed because `cov_mat_b`'s diagonal must include it).
#' @param kc `1` or `3`.
alcd_cv_config <- function(ICC_x, kc = c(1, 3)) {
  kc <- match.arg(as.character(kc), c("1", "3"))

  if (kc == "1") {
    list(
      gamma     = 0.3,
      m_C       = 0,
      ICC_C     = 0.4,
      cov_mat   = matrix(c(1, 0.3,
                            0.3, 1), nrow = 2, byrow = TRUE),
      cov_mat_b = diag(c(ICC_x, 0.4)) + matrix(c(0,   0.1,
                                                  0.1, 0), nrow = 2, byrow = TRUE)
    )
  } else {
    list(
      gamma     = c(0.3, 0.4, 0.25),
      m_C       = c(0, 0, 1),
      ICC_C     = c(0.4, 0.6, 0.35),
      cov_mat   = matrix(c(1,    0.3,  0.2,  0.15,
                            0.3,  1,   -0.1, -0.25,
                            0.2, -0.1,  1,    0.25,
                            0.15, -0.25, 0.25, 1), nrow = 4, byrow = TRUE),
      cov_mat_b = diag(c(ICC_x, 0.4, 0.6, 0.35)) +
        matrix(c(0,    0.1,  0.05,  0.08,
                 0.1,  0,   -0.1,   0.06,
                 0.05, -0.1, 0,    -0.05,
                 0.08,  0.06, -0.05, 0), nrow = 4, byrow = TRUE)
    )
  }
}
