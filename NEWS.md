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