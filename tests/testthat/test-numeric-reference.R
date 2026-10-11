# Regression tests: fixed output on seeded reference datasets. The values were
# generated with dev/make_reference_values.R, which also checks that the
# grid-search optimum of each dataset is well separated from the runner-up
# (relative MSE gap 9.7e-3 for kc = 0 and 2.9e-6 for kc = 3), so they do not
# depend on floating-point details of the platform. If an intended change to
# the estimators alters these values, regenerate them with that script.

test_that("reference output without covariates is unchanged", {
  fit <- reference_fit(kc = 0)

  expect_equal(fit$Coefficients$beta_b, 0.0317119269, tolerance = 1e-7)
  expect_equal(fit$Standard_Error$beta_b, 0.1547118284, tolerance = 1e-7)
  expect_equal(fit$Confidence_Interval$Lower, -0.2715176847, tolerance = 1e-7)
  expect_equal(fit$Confidence_Interval$Upper, 0.3349415384, tolerance = 1e-7)
  expect_equal(fit$Coefficients_ML$beta_b_ML, 0.1271999215, tolerance = 1e-7)
  expect_equal(fit$Standard_Error_ML$beta_b_ML, 1.5523588495, tolerance = 1e-7)
})

test_that("reference output with three covariates is unchanged", {
  fit <- reference_fit(kc = 3)

  expect_equal(unname(unlist(fit$Coefficients)),
               c(0.5809013831, 0.1903223692, 0.2822803030, 0.2207154073), tolerance = 1e-7)
  expect_equal(unname(unlist(fit$Standard_Error)),
               c(0.3145732198, 0.0795971621, 0.0646701376, 0.0794111120), tolerance = 1e-7)
  expect_equal(fit$Confidence_Interval$Lower,
               c(-0.0356507982, 0.0343147982, 0.1555291625, 0.0650724878), tolerance = 1e-7)
  expect_equal(fit$Confidence_Interval$Upper,
               c(1.1974535643, 0.3463299402, 0.4090314435, 0.3763583268), tolerance = 1e-7)
  expect_equal(unname(unlist(fit$Coefficients_ML)),
               c(0.7408268398, 0.1903223692, 0.2822803030, 0.2207154073), tolerance = 1e-7)
  expect_equal(unname(unlist(fit$Standard_Error_ML)),
               c(1.2496096730, 0.0795971621, 0.0646701376, 0.0794111120), tolerance = 1e-7)
})
