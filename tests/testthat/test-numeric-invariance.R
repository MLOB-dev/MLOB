# Exact invariance and equivariance properties of the estimators. These hold
# for any dataset, so they are checked to floating-point precision.

test_that("results do not depend on the row order of the data (no covariates)", {
  df <- reference_data(kc = 0)
  set.seed(1)
  shuffled <- df[sample(nrow(df)), ]
  fit <- mlob(y ~ x, data = shuffled, group = "group")
  expect_equal(fit_values(fit), fit_values(reference_fit(kc = 0)), tolerance = 1e-10)
})

test_that("results do not depend on the row order of the data (kc = 3)", {
  df <- reference_data(kc = 3)
  set.seed(1)
  shuffled <- df[sample(nrow(df)), ]
  fit <- mlob(reference_formula(3), data = shuffled, group = "group")
  expect_equal(fit_values(fit), fit_values(reference_fit(kc = 3)), tolerance = 1e-10)
})

test_that("results do not depend on how groups are labelled", {
  df <- reference_data(kc = 0)
  set.seed(1)
  relabel <- sample(20)
  df$group <- factor(paste0("school_", relabel[df$group]))
  fit <- mlob(y ~ x, data = df, group = "group")
  expect_equal(fit_values(fit), fit_values(reference_fit(kc = 0)), tolerance = 1e-10)
})

test_that("estimates and SE of beta_b are invariant to shifting x and y", {
  df <- reference_data(kc = 0)
  df$x <- df$x + 5
  df$y <- df$y - 100
  fit <- mlob(y ~ x, data = df, group = "group")
  expect_equal(fit_values(fit), fit_values(reference_fit(kc = 0)), tolerance = 1e-10)
})

test_that("rescaling y rescales all estimates and standard errors", {
  df <- reference_data(kc = 3)
  df$y <- 10 * df$y
  fit <- mlob(reference_formula(3), data = df, group = "group")
  expect_equal(fit_values(fit), 10 * fit_values(reference_fit(kc = 3)), tolerance = 1e-10)
})

test_that("rescaling x rescales beta_b and its SE inversely and leaves gamma unchanged", {
  for (kc in c(0, 3)) {
    ref <- reference_fit(kc = kc)
    for (s in c(0.01, 100)) {
      df <- reference_data(kc = kc)
      df$x <- s * df$x
      fit <- mlob(reference_formula(kc), data = df, group = "group")
      expected <- fit_values(ref)
      is_beta <- grepl("^beta_b", names(expected))
      expected[is_beta] <- expected[is_beta] / s
      expect_equal(fit_values(fit), expected, tolerance = 1e-10)
    }
  }
})

test_that("rescaling or shifting a covariate only affects its own coefficient", {
  ref <- reference_fit(kc = 3)
  df <- reference_data(kc = 3)
  df$C1 <- 3 * df$C1
  df$C2 <- df$C2 + 7
  fit <- mlob(reference_formula(3), data = df, group = "group")

  expected <- unlist(ref$Coefficients)
  expected[["gamma_C1"]] <- expected[["gamma_C1"]] / 3
  expect_equal(unlist(fit$Coefficients), expected, tolerance = 1e-10)
  expect_equal(fit$Coefficients_ML$beta_b_ML, ref$Coefficients_ML$beta_b_ML, tolerance = 1e-10)
  expected_se <- unlist(ref$Standard_Error)
  expected_se[["gamma_C1"]] <- expected_se[["gamma_C1"]] / 3
  expect_equal(unlist(fit$Standard_Error), expected_se, tolerance = 1e-10)
  expect_equal(fit$Standard_Error_ML$beta_b_ML, ref$Standard_Error_ML$beta_b_ML, tolerance = 1e-10)
})
