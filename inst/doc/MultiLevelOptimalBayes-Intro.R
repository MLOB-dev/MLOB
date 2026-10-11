## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message    = FALSE,
  warning    = FALSE
)

## ----load-package, echo=TRUE--------------------------------------------------
library(MultiLevelOptimalBayes)

## ----eval = FALSE-------------------------------------------------------------
# mlob(
#   formula,
#   data,
#   group,
#   balancing.limit = 0.2,
#   conf.level = 0.95,
#   jackknife = FALSE,
#   punish.coeff = 2,
#   ...
# )

## -----------------------------------------------------------------------------
result_iris <- mlob(
  Sepal.Length ~ Sepal.Width + Petal.Length,
  data = iris,
  group = "Species",
  conf.level = 0.99,
  jackknife = FALSE
)

summary(result_iris)

## -----------------------------------------------------------------------------
result_chick <- mlob(
  weight ~ Time,
  data = ChickWeight,
  group = "Diet",
  punish.coeff = 1.5,
  jackknife = FALSE
)

print(result_chick)
summary(result_chick)

## -----------------------------------------------------------------------------
result_mtcars <- mlob(
  mpg ~ hp + wt + am + hp:wt + hp:am,
  data = mtcars,
  group = "cyl",
  balancing.limit = 0.35
)

summary(result_mtcars)

## -----------------------------------------------------------------------------
# Get a basic result for demonstration
result <- mlob(weight ~ Time, data = ChickWeight, group = 'Diet', jackknife = FALSE)

# Print method - displays coefficients, standard errors, confidence intervals, Z-values, and p-values
print(result)

## -----------------------------------------------------------------------------
# Summary method - comprehensive summary with significance stars and comparison to unoptimized ML
summary(result)

## -----------------------------------------------------------------------------
# Extract coefficients as a data frame
coef(result)

# Extract standard errors
se(result)

# Extract variance-covariance matrix (diagonal only)
vcov(result)

# Extract confidence intervals
confint(result)

# Extract confidence intervals for specific parameters
confint(result, "beta_b")

# Extract confidence intervals with different confidence level
confint(result, level = 0.99)

## -----------------------------------------------------------------------------
# Convert results to a data frame format
as.data.frame(result)

# Get dimensions (number of parameters)
dim(result)

# Get number of parameters
length(result)

# Get parameter names
names(result)

## -----------------------------------------------------------------------------
# Update model with new parameters (e.g., different confidence level)
updated_result <- update(result, conf.level = 0.99)
summary(updated_result)

## -----------------------------------------------------------------------------
methods(class = "mlob_result")

## -----------------------------------------------------------------------------
# Run analysis
result <- mlob(weight ~ Time, data = ChickWeight, group = 'Diet', jackknife = FALSE)

# Get basic information
cat("Number of parameters:", length(result), "\n")
cat("Parameter names:", paste(names(result), collapse = ", "), "\n")

# Extract key statistics
coefficients <- coef(result)
standard_errors <- se(result)
confidence_intervals <- confint(result, level = 0.99)

# Create a custom summary table
custom_summary <- data.frame(
  Parameter = names(result),
  Estimate = as.numeric(coefficients),
  SE = as.numeric(standard_errors),
  CI_Lower = confidence_intervals[, 1],
  CI_Upper = confidence_intervals[, 2]
)

print(custom_summary)

