# Statistical properties of the regularized Bayes estimator relative to ML,
# as established in Dashuk et al. (2025a): lower standard error, and lower
# mean squared error in small samples with low ICC.

test_that("Bayes SE of beta_b is not larger than the ML SE on the reference datasets", {
  for (kc in c(0, 1, 3)) {
    fit <- reference_fit(kc = kc)
    expect_lte(fit$Standard_Error$beta_b, fit$Standard_Error_ML$beta_b_ML)
  }
})

test_that("Bayes estimator has lower MSE than ML over replications", {
  skip_on_cran()

  beta_b <- 0.5
  est <- t(sapply(101:110, function(seed) {
    set.seed(seed)
    df <- simulate_alcd_cv(J = 20, n = 10, ICC_x = 0.1, beta_b = beta_b, beta_w = 0.2)
    fit <- mlob(y ~ x, data = df, group = "group")
    c(bay = fit$Coefficients$beta_b, ml = fit$Coefficients_ML$beta_b_ML,
      se_bay = fit$Standard_Error$beta_b, se_ml = fit$Standard_Error_ML$beta_b_ML)
  }))

  # Observed on these seeds: RMSE 0.79 (Bayes) vs 1.10 (ML).
  expect_lt(mean((est[, "bay"] - beta_b)^2), mean((est[, "ml"] - beta_b)^2))
  expect_true(all(est[, "se_bay"] <= est[, "se_ml"]))
})
