MultiLevelOptimalBayes
================

# MultiLevelOptimalBayes (MLOB)

**MLOB** is an R package for estimating between-group effects in
multilevel latent variable models using an optimally regularized
Bayesian estimator. It is especially useful for small-sample settings,
low ICC data, and hierarchical models commonly used in psychology,
education, and social sciences.

## Features

- **Regularized Bayesian estimator** optimized for minimum mean squared
  error (MSE)
- **Robust to small sample sizes** and **low intraclass correlation
  (ICC)**
- Supports **covariate adjustment** and **group-level balancing**
- Optional **jackknife resampling** for improved CI coverage
- Returns **full inferential output**: point estimates, standard errors,
  p-values, and confidence intervals
- Clean and user-friendly interface via the `mlob()` function
- **Comprehensive S3 methods** following standard R conventions

## Installation

To install the development version from GitHub:

``` r
install.packages("devtools")
devtools::install_github("MLOB-dev/MLOB")
```

MLOB is available on CRAN under the GPL-3 license. To install the
released version:

``` r
install.packages("MultiLevelOptimalBayes")
```

## View the Vignette

After installing the package, run the following to open the introductory
vignette:

``` r
vignette("MultiLevelOptimalBayes-Intro")
```

## Examples

``` r
library(MultiLevelOptimalBayes)
```

### Fit a model on the iris dataset

``` r
result <- mlob(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris,
               group = "Species", conf.level = 0.95)
```

### View results

``` r
summary(result)
#> Call:
#>  mlob(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris, group = Species, conf.level = 0.95) 
#> 
#> Summary of Coefficients:
#>                     Estimate Std. Error Lower CI (95%) Upper CI (95%)   Z value
#> beta_b             0.8308711  1.4655556    -2.04156502      3.7033072 0.5669325
#> gamma_Petal.Length 0.4679522  0.2582579    -0.03822406      0.9741285 1.8119567
#>                      Pr(>|z|) Significance
#> beta_b             0.57076004             
#> gamma_Petal.Length 0.06999289            .
#> 
#> 
#> For comparison, summary of coefficients from unoptimized analysis (ML):
#>                     Estimate   Std. Error Lower CI (95%) Upper CI (95%)
#> beta_b             0.6027440 5.424780e+15  -1.063237e+16   1.063237e+16
#> gamma_Petal.Length 0.4679522 2.582579e-01  -3.822406e-02   9.741285e-01
#>                         Z value   Pr(>|z|) Significance
#> beta_b             1.111094e-16 1.00000000             
#> gamma_Petal.Length 1.811957e+00 0.06999289            .
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Note:
#>   The standard error from unoptimized ML estimation is about 3.701518e+17% larger than the standard error obtained through our optimization procedure,
#>   meaning that the optimized estimates are more accurate.
#>   Concerning the estimates themselves, the unoptimized ML estimates may
#>   differ greatly from the optimized estimates and should not be reported.
#>   As the optimized estimates are always at least as accurate as the
#>   unoptimized ML estimates,
#>   please use them and their corresponding standard errors (first table of
#>   output) for interpretation and reporting.
#>   For more information, see Dashuk et al. (2025).
```

## Available Methods

The `mlob_result` object supports a comprehensive set of methods that
follow standard R conventions:

### Display Methods

``` r
print(result)      # Display coefficients, standard errors, confidence intervals, Z-values, and p-values
#> Call:
#>  mlob(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris, group = Species, conf.level = 0.95) 
#> 
#> Coefficients
#>     beta_b gamma_Petal.Length
#>  0.8308711          0.4679522
#> 
#> Standard_Error
#>    beta_b gamma_Petal.Length
#>  1.465556          0.2582579
#> 
#> Confidence_Interval (95%)
#>                          Lower     Upper
#> beta_b             -2.04156502 3.7033072
#> gamma_Petal.Length -0.03822406 0.9741285
#> 
#> Z value
#>     beta_b gamma_Petal.Length
#>  0.5669325           1.811957
#> 
#> p value
#>   beta_b gamma_Petal.Length
#>  0.57076         0.06999289
```

``` r
summary(result)   # Comprehensive summary with significance stars and comparison to unoptimized ML
#> Call:
#>  mlob(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris, group = Species, conf.level = 0.95) 
#> 
#> Summary of Coefficients:
#>                     Estimate Std. Error Lower CI (95%) Upper CI (95%)   Z value
#> beta_b             0.8308711  1.4655556    -2.04156502      3.7033072 0.5669325
#> gamma_Petal.Length 0.4679522  0.2582579    -0.03822406      0.9741285 1.8119567
#>                      Pr(>|z|) Significance
#> beta_b             0.57076004             
#> gamma_Petal.Length 0.06999289            .
#> 
#> 
#> For comparison, summary of coefficients from unoptimized analysis (ML):
#>                     Estimate   Std. Error Lower CI (95%) Upper CI (95%)
#> beta_b             0.6027440 5.424780e+15  -1.063237e+16   1.063237e+16
#> gamma_Petal.Length 0.4679522 2.582579e-01  -3.822406e-02   9.741285e-01
#>                         Z value   Pr(>|z|) Significance
#> beta_b             1.111094e-16 1.00000000             
#> gamma_Petal.Length 1.811957e+00 0.06999289            .
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Note:
#>   The standard error from unoptimized ML estimation is about 3.701518e+17% larger than the standard error obtained through our optimization procedure,
#>   meaning that the optimized estimates are more accurate.
#>   Concerning the estimates themselves, the unoptimized ML estimates may
#>   differ greatly from the optimized estimates and should not be reported.
#>   As the optimized estimates are always at least as accurate as the
#>   unoptimized ML estimates,
#>   please use them and their corresponding standard errors (first table of
#>   output) for interpretation and reporting.
#>   For more information, see Dashuk et al. (2025).
```

### Statistical Methods

``` r
coef(result)                    # Extract coefficients as a data frame
#>      beta_b gamma_Petal.Length
#> 1 0.8308711          0.4679522
se(result)                      # Extract standard errors
#>             beta_b gamma_Petal.Length 
#>          1.4655556          0.2582579
vcov(result)                    # Extract variance-covariance matrix (diagonal only)
#>             beta_b gamma_Petal.Length 
#>         2.14785309         0.06669717
confint(result)                 # Extract confidence intervals
#>                           2.5%     97.5%
#> beta_b             -2.04156502 3.7033072
#> gamma_Petal.Length -0.03822406 0.9741285
confint(result, "beta_b")       # Extract CI for specific parameter
#>             2.5%    97.5%
#> beta_b -2.041565 3.703307
confint(result, level = 0.99)   # Extract CI with different confidence level
#>                          0.5%    99.5%
#> beta_b             -2.9441499 4.605892
#> gamma_Petal.Length -0.1972762 1.133181
```

### Utility Methods

``` r
as.data.frame(result)           # Convert results to a data frame format
#>                     Estimate Std. Error Lower CI (95%) Upper CI (95%)   Z value
#> beta_b             0.8308711  1.4655556    -2.04156502      3.7033072 0.5669325
#> gamma_Petal.Length 0.4679522  0.2582579    -0.03822406      0.9741285 1.8119567
#>                      Pr(>|z|)
#> beta_b             0.57076004
#> gamma_Petal.Length 0.06999289
dim(result)                     # Get dimensions (number of parameters)
#> [1] 1 2
length(result)                  # Get number of parameters
#> [1] 2
names(result)                   # Get parameter names
#> [1] "beta_b"             "gamma_Petal.Length"
```

### Update Method

``` r
updated_result <- update(result, conf.level = 0.99)  # Update model with new parameters
summary(updated_result)
#> Call:
#>  mlob(Sepal.Length ~ Sepal.Width + Petal.Length, data = data, group = Species, conf.level = 0.99, jackknife = FALSE) 
#> 
#> Summary of Coefficients:
#>                     Estimate Std. Error Lower CI (99%) Upper CI (99%)   Z value
#> beta_b             0.8308711  1.4655556     -2.9441499       4.605892 0.5669325
#> gamma_Petal.Length 0.4679522  0.2582579     -0.1972762       1.133181 1.8119567
#>                      Pr(>|z|) Significance
#> beta_b             0.57076004             
#> gamma_Petal.Length 0.06999289            .
#> 
#> 
#> For comparison, summary of coefficients from unoptimized analysis (ML):
#>                     Estimate   Std. Error Lower CI (99%) Upper CI (99%)
#> beta_b             0.6027440 5.424780e+15  -1.397331e+16   1.397331e+16
#> gamma_Petal.Length 0.4679522 2.582579e-01  -1.972762e-01   1.133181e+00
#>                         Z value   Pr(>|z|) Significance
#> beta_b             1.111094e-16 1.00000000             
#> gamma_Petal.Length 1.811957e+00 0.06999289            .
#> 
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Note:
#>   The standard error from unoptimized ML estimation is about 3.701518e+17% larger than the standard error obtained through our optimization procedure,
#>   meaning that the optimized estimates are more accurate.
#>   Concerning the estimates themselves, the unoptimized ML estimates may
#>   differ greatly from the optimized estimates and should not be reported.
#>   As the optimized estimates are always at least as accurate as the
#>   unoptimized ML estimates,
#>   please use them and their corresponding standard errors (first table of
#>   output) for interpretation and reporting.
#>   For more information, see Dashuk et al. (2025).
```

### Discover Available Methods

``` r
methods(class = "mlob_result")  # List all available methods
#>  [1] as.data.frame coef          confint       dim           length       
#>  [6] names         print         se            summary       update       
#> [11] vcov         
#> see '?methods' for accessing help and source code
```

All methods follow standard R conventions, making `mlob_result` objects
compatible with existing R workflows and familiar to users of other
statistical packages.

## Limitations

- The estimator assumes approximately equal group sizes. Although
  balancing helps, unequal sizes may still bias results.

- Grid-search is local around the ML estimate; global optimum is found
  with high probability but is not guaranteed.

- Jackknife resampling improves inference in small samples but can be
  computationally heavy in larger samples.

- Currently supports two-level models with continuous outcomes only.
  Extensions to GLMMs or 3+ level models are future work.

## Contributing & Support

Please open an issue at:

<https://github.com/MLOB-dev/MLOB/issues>

Users may also join discussions or suggest enhancements on the
Discussions page at

<https://github.com/MLOB-dev/MLOB/discussions>.

## Authors

Valerii Dashuk

Binayak Timilsina

Martin Hecht

Steffen Zitzmann

## Citation

If you use MLOB in your research, please cite:

Dashuk, V., Hecht, M., Luedtke, O., Robitzsch, A., & Zitzmann, S.
(2024). An Optimally Regularized Estimator of Multilevel Latent Variable
Models, with Improved MSE Performance
<https://doi.org/10.13140/RG.2.2.18148.39048>

## Contact:

<martin.hecht@hsu-hh.de>

<steffen.zitzmann@medicalschool-hamburg.de>

<multilob@outlook.com>
