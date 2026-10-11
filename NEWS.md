# MultiLevelOptimalBayes 0.1.0.0

This release corrects the implementation of the regularized Bayesian
estimator so that it follows the derivation in Dashuk et al. (2025a, 2025b).
**Estimates and standard errors of `beta_b` change**; results published with
earlier versions can be reproduced with MultiLevelOptimalBayes 0.0.5.0 or
earlier from the CRAN archive.

## Bug fixes in the Bayesian estimator

* The gamma approximations of the estimated between-group variance of the
  predictor (`tau_x2`) and covariance (`tau_yx`) had the wrong mean: the
  eigenvectors of the within-group part were never filled in, so the
  approximations described the between-group mean square `MSA / n`
  (mean `tau + sigma / n`) instead of `(MSA - MSD) / n` (mean `tau`). The mean
  and variance of both estimators are now computed exactly in closed form.
* The scale parameters of the gamma approximations lacked a factor of 2
  (each `theta_i * chi^2_1` term is `Gamma(1/2, 2 theta_i)`), which halved
  their means and made the "5 standard deviation" search window 2.5 standard
  deviations wide.
* The MSE minimized in the grid search now is the one derived in
  Dashuk et al. (2025a, Eqs. 78 and 101). The previous implementation had the
  mean of the denominator instead of its square in the first term, which made
  the estimate depend on the units of the predictor.
* The standard errors of `beta_b` (Bayesian and ML) use the same corrected
  moments. Where the variance of the ML ratio does not exist under the gamma
  approximation (`k_B <= 2`), a first-order approximation is reported.
* Standard errors of the covariate coefficients (`gamma`) are now the OLS
  standard errors of the regression that estimates them. Previously the
  residuals omitted the intercept and the other covariates, so the standard
  errors depended on the mean of `y` and on the location of the covariates.

## Other changes

* Estimates of `beta_b` are now exactly equivariant to the units of `x` and
  `y` (multiplying `x` by `s` divides `beta_b` and its standard error by `s`).
* Estimation no longer builds `(nJ + J + 1)`-dimensional covariance matrices,
  so `mlob()` (and especially `jackknife = TRUE`) is much faster and works for
  large data sets. The `pracma` dependency is removed.
* Example 1 in `?mlob` used `conf.level = 0.01` (a 1% interval); it now uses
  `0.99`.

# MultiLevelOptimalBayes 0.0.5.0

## Bug fixes

* `mlob()` now gives correct results when the rows of a group are not
  contiguous in `data`. Previously the estimators silently assigned
  observations to the wrong groups unless the data were sorted by the grouping
  variable, which changed the estimates (e.g. the `mtcars` example in `?mlob`,
  grouped by the unsorted `cyl`).
* `mlob()` now warns when the estimated between-group variance of the
  predictor is negative. In that case the estimates of `beta_b` are not
  meaningful; previously they were returned without notice.
* The checks for missing between- and within-group variation now treat
  (co)variance components as zero up to a tolerance relative to the variances
  of the data, instead of testing for exact equality with zero, which could be
  missed because of floating-point rounding.

## Testing

* Added a simulation-based test fixture (`simulate_alcd_cv()`) implementing the
  data-generating process published in the supplementary material of
  Dashuk et al. (2025a, 2025b), for the no-covariate model and any number of
  control variables.
* Added numeric tests on seeded simulated data:
  * ML estimates agree with an independent `anova(lm())` implementation of the
    estimator (with and without control variables) and with variance
    components from REML fits in `lme4`.
  * The Bayes estimate equals the regularized ratio at the grid-search optimum,
    and confidence intervals and p-values are consistent with the reported
    estimates and standard errors.
  * Results are invariant to row order, group labels and shifts of the
    variables, and equivariant to rescaling of the response and covariates.
  * Reference values guard against unintended changes in output.
  * The Bayes estimator has a lower standard error and lower mean squared error
    than ML across replications.

## Documentation

* The `ChickWeight` example in `?mlob`, which demonstrates all S3 methods, is
  now wrapped in `\donttest{}` to keep example run time within CRAN limits.
* Fixed a stray comment marker that caused `#'` to be rendered literally in the
  `?mlob` examples section.

## Internal

* Added `lme4` to `Suggests` (used in tests only).
