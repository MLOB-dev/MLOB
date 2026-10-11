# MultiLevelOptimalBayes 0.1.0.0

This release corrects the implementation of the regularized Bayesian
estimator so that it follows the derivation in Dashuk et al. (2025a, 2025b).
**Estimates and standard errors of `beta_b` change**; results from earlier
versions can be reproduced with MultiLevelOptimalBayes 0.0.4.0 from the CRAN
archive. (GitHub branch `v0.0.5.0` has the same estimator family, with
additional bug fixes.)

## Bug fixes in the Bayesian estimator

* The gamma approximations described the between-group mean square `MSA / n`
  (mean `tau + sigma / n`) instead of `(MSA - MSD) / n` (mean `tau`). The mean
  and variance of both estimators are now computed exactly in closed form.
* Changed the scale parameters of the gamma approximations to improve the
  search window.
* Fixed dependence of the estimate on the units of the predictor.
* The standard errors of `beta_b` (Bayesian and ML) use the same corrected
  moments. Where the variance of the ML ratio does not exist under the gamma
  approximation (`k_B <= 2`), a first-order approximation is reported.
* Improved the computation of the standard errors of the covariate
  coefficients (`gamma`).

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
  contiguous in `data`.
* `mlob()` now gives an explicit warning when the estimated between-group
  variance of the predictor is negative. Previously this was only hinted at
  later by messages such as "Please check grouping" or "Singularities in the
  between-group estimation".
* Added an early check for near-zero between- and within-group (co)variances.

## Testing

* Added a simulation-based test fixture (`simulate_alcd_cv()`) for the test
  suite only (not exported; implements the data-generating process from the
  supplementary material of Dashuk et al. (2025a, 2025b)).
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
