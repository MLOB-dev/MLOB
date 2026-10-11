# Seeded reference datasets shared by the numeric test files. Fits are cached
# so each reference model is estimated only once per test run.
#
# Settings: J = 20 groups of n = 10, ICC_x = 0.1, beta_w = 0.2, ICC_y = 0.2;
# beta_b = 0.5 without control variables, beta_b = 0.2 for kc = 1 and kc = 3
# (with the ALCD_CV.m covariance constants, larger beta_b gives a negative
# residual variance for kc = 3).

reference_data <- function(kc = 0, seed = 2) {
  set.seed(seed)
  if (kc == 0) {
    simulate_alcd_cv(J = 20, n = 10, ICC_x = 0.1, beta_b = 0.5, beta_w = 0.2)
  } else {
    do.call(simulate_alcd_cv,
            c(list(J = 20, n = 10, ICC_x = 0.1, beta_b = 0.2, beta_w = 0.2),
              alcd_cv_config(0.1, kc)))
  }
}

reference_formula <- function(kc = 0) {
  if (kc == 0) return(y ~ x)
  stats::as.formula(paste("y ~ x +", paste0("C", seq_len(kc), collapse = " + ")))
}

.reference_cache <- new.env()

reference_fit <- function(kc = 0, seed = 2) {
  key <- paste(kc, seed)
  if (is.null(.reference_cache[[key]])) {
    .reference_cache[[key]] <- mlob(reference_formula(kc), data = reference_data(kc, seed),
                                    group = "group")
  }
  .reference_cache[[key]]
}

# All numeric output of an mlob_result that should be invariant (or
# equivariant) under data transformations.
fit_values <- function(fit) {
  c(unlist(fit$Coefficients), unlist(fit$Standard_Error),
    unlist(fit$Coefficients_ML), unlist(fit$Standard_Error_ML))
}
