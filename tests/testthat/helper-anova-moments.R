# Independent reference implementation of the latent-covariate ML estimator
# (Luedtke et al. 2008), built on stats::anova(lm()) instead of the package's
# hand-written sums of squares. Cross-products are obtained by polarization:
# cov(u, v) = (var(u + v) - var(u) - var(v)) / 2, applied to the one-way
# ANOVA variance components.

anova_components <- function(v, g) {
  ms <- stats::anova(stats::lm(v ~ factor(g)))[["Mean Sq"]]
  n <- length(v) / length(unique(g))
  c(tau = (ms[1] - ms[2]) / n, sigma = ms[2])
}

anova_cross_components <- function(u, v, g) {
  (anova_components(u + v, g) - anova_components(u, g) - anova_components(v, g)) / 2
}

#' ML estimate of beta_b (and gamma) from `simulate_alcd_cv()`-style data,
#' following the estimator of Dashuk et al. (2025b) for control variables.
anova_ml_estimate <- function(df, kc = 0) {
  g <- df$group
  vx <- anova_components(df$x, g)
  y <- df$y
  gamma <- se_gamma <- NULL
  if (kc > 0) {
    C <- as.matrix(df[paste0("C", seq_len(kc))])
    phi <- sapply(seq_len(kc), function(i) anova_cross_components(df$x, C[, i], g) / vx)
    x_b <- stats::ave(df$x, g)
    C_tilde <- C - cbind(x_b, df$x - x_b) %*% phi
    coefs <- summary(stats::lm(y ~ C_tilde))$coefficients
    gamma <- unname(coefs[-1, "Estimate"])
    se_gamma <- unname(coefs[-1, "Std. Error"])
    y <- as.numeric(y - C %*% gamma)
  }
  vxy <- anova_cross_components(df$x, y, g)
  list(beta_b_ML = unname(vxy[["tau"]] / vx[["tau"]]), gamma = gamma, se_gamma = se_gamma,
       tau_x2 = unname(vx[["tau"]]), tau_yx = unname(vxy[["tau"]]))
}

#' Covariance matrix of z = (v_11, ..., v_nJ, group means, grand mean) under the
#' random-intercept model with between component `tau` and within component
#' `sigma` (Dashuk et al. 2025a, Appendix F). The ANOVA estimator of the
#' between component is the quadratic form z' A z with `A = anova_weights(n, J)`.
latent_cov <- function(tau, sigma, n, J) {
  N <- n * J; G <- rep(seq_len(J), each = n); p <- N + J + 1
  S <- matrix(0, p, p)
  S[1:N, 1:N] <- outer(G, G, "==") * tau + diag(sigma, N)
  S[1:N, N + seq_len(J)] <- outer(G, seq_len(J), "==") * (tau + sigma / n)
  S[1:N, p] <- tau / J + sigma / (n * J)
  S[N + seq_len(J), N + seq_len(J)] <- diag(tau + sigma / n, J)
  S[N + seq_len(J), p] <- tau / J + sigma / (n * J)
  S[p, p] <- tau / J + sigma / (n * J)
  S[lower.tri(S)] <- t(S)[lower.tri(S)]
  S
}

anova_weights <- function(n, J) {
  diag(c(rep(-1 / (n * (n - 1) * J), n * J),
         rep((n * J - 1) / ((n - 1) * (J - 1) * J), J),
         -J / (J - 1)))
}
