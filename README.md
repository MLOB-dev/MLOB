
# MultiLevelOptimalBayes

<!-- badges: start -->
<!-- badges: end -->

# MultiLevelOptimalBayes (MLOB) <img src="https://img.shields.io/badge/R-CRAN-blue.svg" alt="R badge" />

**MLOB** is an R package for estimating between-group effects in
multilevel latent variable models using an optimally regularized
Bayesian estimator. It is especially useful for small-sample settings,
low ICC data, and hierarchical models commonly used in psychology,
education, and social sciences.

------------------------------------------------------------------------

## ✨ Features

- ✅ **Regularized Bayesian estimator** optimized for minimum mean
  squared error (MSE)
- ✅ **Robust to small sample sizes** and **low intraclass correlation
  (ICC)**
- ✅ Supports **covariate adjustment** and **group-level balancing**
- ✅ Optional **jackknife resampling** for improved CI coverage
- ✅ Returns **full inferential output**: point estimates, standard
  errors, p-values, and confidence intervals
- ✅ Clean and user-friendly interface via the `mlob()` function

------------------------------------------------------------------------

## 📦 Installation

To install the development version from GitHub:

\`\`\`r install.packages(“devtools”)

MLOB is available on CRAN under the GPL-3 license. To install the
released version: \`\`\`r install.packages(“MLOB”)

The development version can be installed from GitHub via devtools:
\`\`\`r devtools::install_github(“Dfkthbq01/MultiLevelOptimalBayes”)

## 📦 View the Vignette

After installing the package, run the following to open the introductory
vignette:

\`\`\`r vignette(“MultiLevelOptimalBayes-Intro”)

## 📦 Examples

\`\`\`r

library(MultiLevelOptimalBayes)

### Fit a model on the iris dataset

result \<- mlob( Sepal.Length ~ Sepal.Width + Petal.Length, data = iris,
group = “Species”, conf.level = 0.95 )

### View results

summary(result)

## 📦 Limitations

-The estimator assumes approximately equal group sizes. Although
balancing helps, unequal sizes may still bias results.

- Grid-search is local (±5σ) around the ML estimate; global optimum is
  not guaranteed.

- Linear covariate residualization may miss nonlinear effects or
  interactions.

- Jackknife resampling improves inference in small samples but can be
  computationally heavy.

-Currently supports two-level models with continuous outcomes only.
Extensions to GLMMs or 3+ level models are future work.

## 📦 Contributing & Support

Please open an issue at:

<https://github.com/Dfkthbq01/MultiLevelOptimalBayes/issues>

Users may also join discussions or suggest enhancements on the
Discussions page at

<https://github.com/Dfkthbq01/MultiLevelOptimalBayes/discussions>.

## 📦Authors

Valerii Dashuk

Binayak Timilsina

Martin Hecht

Steffen Zitzmann

## 📚 Citation

If you use MLOB in your research, please cite:

Dashuk, V., Hecht, M., Luedtke, O., Robitzsch, A., & Zitzmann, S.
(2024). An Optimally Regularized Estimator of Multilevel Latent Variable
Models with Improved MSE Performance.
<https://doi.org/10.13140/RG.2.2.18148.39048>

## 📫 Contact:

<martin.hecht@hsu-hh.de>

<steffen.zitzmann@medicalschool-hamburg.de>

<multilob@outlook.com>
