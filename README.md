# PSPI
  <!-- badges: start -->
  [![R-CMD-check](https://github.com/zjg540066169/SBMtrees/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zjg540066169/PSPI/actions/workflows/R-CMD-check.yaml)
  [![License: GPL-2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
  [![version](https://img.shields.io/badge/version-1.1-green.svg)](https://github.com/zjg540066169/PSPI)
  ![R](https://img.shields.io/badge/language-R-blue)
  ![C++](https://img.shields.io/badge/language-C%2B%2B-green)
  <!-- badges: end -->

The R package **PSPI** (Propensity Score Predictive Inference) provides Bayesian methods for generalizing treatment effects from clinical trials to target populations. It implements five models-\code{BCF}, \code{BCF_P}, \code{FullBART}, \code{SplineBART}, and \code{DSplineBART}-built on Bayesian
Additive Regression Trees (BART). Spline-based variants (\code{SplineBART} and \code{DSplineBART}) use propensity score transformations and spline terms to handle covariate shift between datasets.
Core computations rely on efficient MCMC routines implemented in C++.

## Installation
This package is based on `Rcpp`, `RcppArmadillo`, `RcppDist`, and `pg`, please make sure these three packages can be installed.

This package can be installed from R CRAN:
```
install.packages("PSPI")
```
or Github:
```
require("devtools")
install_github("https://github.com/zjg540066169/PSPI")
library(PSPI)
```

## Attribution

This package includes code derived from the [BART3](https://github.com/rsparapa/bnptools/tree/master) package, originally developed by Rodney Sparapani. 

The original source code, licensed under the [GNU General Public License version 2 (GPL-2)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html), has been modified as follows:
- We include part of the C++ code in BART3, primarily about functions about `wbart` and `cpwart`. We also modify some files to make sure our package can be successfully compiled.
- Modifications were made by Jungang Zou, 2024.

### Licensing

- The original BART3 package is licensed under the GNU General Public License version 2 (GPL-2).
- This package, as a derived work, is also licensed under the GNU General Public License version 2 (GPL-2) to comply with the licensing terms.



