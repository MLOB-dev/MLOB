# Numeric correctness of the estimators against independent implementations.
# Data: seeded draws from the published ALCD / ALCD_CV data-generating process
# (see helper-simulate-alcd.R and helper-reference-fits.R).

test_that("ML estimate equals an independent anova(lm()) computation (no covariates)", {
  df <- reference_data(kc = 0)
  fit <- reference_fit(kc = 0)
  ref <- anova_ml_estimate(df)

  expect_equal(fit$Coefficients_ML$beta_b_ML, ref$beta_b_ML, tolerance = 1e-10)
})

test_that("ML estimates equal an independent anova(lm()) computation (kc = 1 and kc = 3)", {
  for (kc in c(1, 3)) {
    df <- reference_data(kc = kc)
    fit <- reference_fit(kc = kc)
    ref <- anova_ml_estimate(df, kc = kc)

    expect_equal(fit$Coefficients_ML$beta_b_ML, ref$beta_b_ML, tolerance = 1e-10)
    expect_equal(unname(unlist(fit$Coefficients_ML[-1])), ref$gamma, tolerance = 1e-10)
    expect_equal(unname(unlist(fit$Coefficients[-1])), ref$gamma, tolerance = 1e-10)
    expect_equal(unname(unlist(fit$Standard_Error[-1])), ref$se_gamma, tolerance = 1e-10)
  }
})

test_that("moments of tau_x2 and tau_yx equal the exact moments of the quadratic forms", {
  # tau_hat = z' A z with z ~ N(0, Sigma): E = tr(A Sigma), Var = 2 tr((A Sigma)^2); for
  # the cross product z_x' A z_y: Var = tr(A S_x A S_y) + tr((A S_xy)^2).
  for (seed in 2:4) {
    df <- reference_data(kc = 0, seed = seed)
    d <- list(y = df$y, x = df$x, k = 20, n = 10, kn = 200)
    d$b_b <- estimate_ML_CV(d)$beta_b_ML_CV
    bay <- estimate_Bay_CV(d)
    A <- anova_weights(10, 20)
    Mx <- A %*% latent_cov(bay$tau_x2, bay$sigma_x2, 10, 20)
    My <- A %*% latent_cov(bay$tau_y2, bay$sigma_y2, 10, 20)
    Mxy <- A %*% latent_cov(bay$tau_yx, bay$sigma_yx, 10, 20)

    expect_equal(bay$E_x, sum(diag(Mx)), tolerance = 1e-10)
    expect_equal(bay$V_x, 2 * sum(Mx * t(Mx)), tolerance = 1e-10)
    expect_equal(bay$E_yx, sum(diag(Mxy)), tolerance = 1e-10)
    expect_equal(bay$V_yx, sum(Mx * t(My)) + sum(Mxy * t(Mxy)), tolerance = 1e-10)
  }
})

test_that("ML estimate equals the lme4 REML variance-component reconstruction", {
  skip_if_not_installed("lme4")

  # For a balanced one-way design, REML variance components equal the ANOVA
  # estimates whenever those are positive, so beta_b_ML = tau_xy / tau_x2 can
  # be rebuilt from random-intercept fits of x, y and x + y.
  for (seed in c(2, 3)) {
    df <- reference_data(kc = 0, seed = seed)
    tau_lmer <- function(v) {
      m <- lme4::lmer(v ~ 1 + (1 | g), data = data.frame(v = v, g = df$group), REML = TRUE)
      expect_false(lme4::isSingular(m))
      as.numeric(lme4::VarCorr(m)$g)
    }
    tau_x <- tau_lmer(df$x)
    tau_xy <- (tau_lmer(df$x + df$y) - tau_x - tau_lmer(df$y)) / 2

    expect_equal(reference_fit(kc = 0, seed = seed)$Coefficients_ML$beta_b_ML,
                 tau_xy / tau_x, tolerance = 1e-5)
  }
})

test_that("Bayes estimate is the regularized ratio at the grid-search optimum", {
  for (kc in c(0, 3)) {
    df <- reference_data(kc = kc)
    fit <- reference_fit(kc = kc)
    d <- list(y = df$y, x = df$x, k = 20, n = 10, kn = 200)
    if (kc > 0) {
      d$C <- as.matrix(df[paste0("C", seq_len(kc))])
      d$kc <- kc
    }
    d$b_b <- estimate_ML_CV(d)$beta_b_ML_CV
    bay <- estimate_Bay_CV(d)

    opt <- which.min(bay$MSE)
    w <- bay$W[opt]
    tau02 <- bay$Tau02[opt]
    expect_true(w > 0 && w <= 1)
    expect_true(tau02 >= 0.01)
    expect_equal(bay$beta_b_Bay, bay$tau_yx / ((1 - w) * tau02 + w * bay$tau_x2),
                 tolerance = 1e-12)
    expect_equal(fit$Coefficients$beta_b, bay$beta_b_Bay, tolerance = 1e-12)
    expect_equal(fit$Coefficients_ML$beta_b_ML, bay$tau_yx / bay$tau_x2, tolerance = 1e-12)

    ref <- anova_ml_estimate(df, kc = kc)
    expect_equal(bay$tau_x2, ref$tau_x2, tolerance = 1e-10)
    expect_equal(bay$tau_yx, ref$tau_yx, tolerance = 1e-10)
  }
})

test_that("confidence intervals are estimate +/- qnorm * SE at the requested level", {
  fit <- reference_fit(kc = 3)
  z <- stats::qnorm(0.975)
  est <- unlist(fit$Coefficients)
  se <- unlist(fit$Standard_Error)
  expect_equal(fit$Confidence_Interval$Lower, unname(est - z * se), tolerance = 1e-12)
  expect_equal(fit$Confidence_Interval$Upper, unname(est + z * se), tolerance = 1e-12)

  fit99 <- mlob(reference_formula(3), data = reference_data(kc = 3), group = "group",
                conf.level = 0.99)
  z99 <- stats::qnorm(0.995)
  expect_equal(unlist(fit99$Coefficients), est, tolerance = 1e-12)
  expect_equal(fit99$Confidence_Interval$Upper - fit99$Confidence_Interval$Lower,
               unname(2 * z99 * se), tolerance = 1e-12)
  expect_equal(unlist(fit$p_value), 2 * stats::pnorm(-abs(est / se)), tolerance = 1e-10)
})
