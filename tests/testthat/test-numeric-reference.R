# Regression tests: fixed output on seeded reference datasets. The values were
# generated with dev/make_reference_values.R, which also checks that the
# grid-search optimum of each dataset is well separated from the runner-up
# (relative MSE gap 2.5e-3 for kc = 0 and 4.7e-4 for kc = 3), so they do not
# depend on floating-point details of the platform. If an intended change to
# the estimators alters these values, regenerate them with that script.

test_that("reference output without covariates is unchanged", {
  fit <- reference_fit(kc = 0)

  expect_equal(fit$Coefficients$beta_b, 0.0511055186, tolerance = 1e-7)
  expect_equal(fit$Standard_Error$beta_b, 0.1242312232, tolerance = 1e-7)
  expect_equal(fit$Confidence_Interval$Lower, -0.1923832047, tolerance = 1e-7)
  expect_equal(fit$Confidence_Interval$Upper, 0.2945942419, tolerance = 1e-7)
  expect_equal(fit$Coefficients_ML$beta_b_ML, 0.1271999215, tolerance = 1e-7)
  expect_equal(fit$Standard_Error_ML$beta_b_ML, 0.4014984292, tolerance = 1e-7)
})

test_that("reference output with three covariates is unchanged", {
  fit <- reference_fit(kc = 3)

  expect_equal(unname(unlist(fit$Coefficients)),
               c(0.3129944773, 0.1903223692, 0.2822803030, 0.2207154073), tolerance = 1e-7)
  expect_equal(unname(unlist(fit$Standard_Error)),
               c(0.0843644834, 0.0799086265, 0.0692037295, 0.0746493502), tolerance = 1e-7)
  expect_equal(fit$Confidence_Interval$Lower,
               c(0.1476431284, 0.0337043392, 0.1466434855, 0.0744053695), tolerance = 1e-7)
  expect_equal(fit$Confidence_Interval$Upper,
               c(0.4783458263, 0.3469403991, 0.4179171204, 0.3670254451), tolerance = 1e-7)
  expect_equal(unname(unlist(fit$Coefficients_ML)),
               c(0.7408268398, 0.1903223692, 0.2822803030, 0.2207154073), tolerance = 1e-7)
  expect_equal(unname(unlist(fit$Standard_Error_ML)),
               c(0.3490383924, 0.0799086265, 0.0692037295, 0.0746493502), tolerance = 1e-7)
})
