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
  gamma <- NULL
  if (kc > 0) {
    C <- as.matrix(df[paste0("C", seq_len(kc))])
    phi <- sapply(seq_len(kc), function(i) anova_cross_components(df$x, C[, i], g) / vx)
    x_b <- stats::ave(df$x, g)
    C_tilde <- C - cbind(x_b, df$x - x_b) %*% phi
    gamma <- unname(stats::coef(stats::lm(y ~ C_tilde))[-1])
    y <- as.numeric(y - C %*% gamma)
  }
  vxy <- anova_cross_components(df$x, y, g)
  list(beta_b_ML = unname(vxy[["tau"]] / vx[["tau"]]), gamma = gamma,
       tau_x2 = unname(vx[["tau"]]), tau_yx = unname(vxy[["tau"]]))
}
