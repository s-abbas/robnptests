
robTests
========

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/s-abbas/robTests/branch/master/graph/badge.svg)](https://codecov.io/gh/s-abbas/robTests?branch=master)
<!-- badges: end -->

The `robTests` R-package contains different robust and non-parametric tests for the two-sample location problem. The tests allow for comparisons of either the means or the scales of two random samples.

Installation
------------

To install the package, the `devtools` package is required.

``` r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github("s-abbas/robTests")

library(robTests)
```

Scope and Usage
---------------

The robust and non-parametric tests contained in this R-package follow the construction principle of the popular t-test: An estimate for the location difference of the two samples is divided by an estimate of scale. Both estimators are robust. The p-values can either be computed using the permutation principle, the randomization principle, or the asymptotic distribution of the estimators. An appropriate test principle is selected automatically or may be specified by the user. 
Functions that implement the estimators of scale and location used by the tests are made available as well.

The implemented tests are:

* tests based on the median (`med_test`) and one- and two-sample Hodges-Lehmann estimators (`hl1_test()`, resp. `hl2_test()`), scaled by appropriate robust estimators
* tests based on the trimmed mean and winsorized variance (`trimmed_test()`) 
* tests based on the Huber-, Hampel- or Bisquare-M-estimator (`m_test()`), if necessary divided by a pooled tau-estimate

Tests for a difference in scale are performed by transforming the two samples prior to applying the test. For details on the tests and references see the documentation of the different functions as well as the accompanying vignettes.

### Example: Performing the `hl2_test`

``` r
set.seed(121)
x <- rnorm(50); y <- rnorm(50)

## Perform an asymptotic test based on the two sample Hodges-Lehmann estimator
hl2_test(x, y, method = "asymptotic")

# 
# 	Asymptotic test based on HL2-estimator
# 
# data:  x and y
# D = 1.0916, p-value = 0.275
# alternative hypothesis: true location shift is not equal to 0
# sample estimates:
# HL2 of x and y 
#      0.2048249

## Perform an according asymptotic test for a difference in scale
hl2_test(x, y, method = "asymptotic", var.test = TRUE)

# 	Asymptotic test based on HL2-estimator
# 
# data:  x and y
# S = -0.24094, p-value = 0.8096
# alternative hypothesis: true ratio of variances is not equal to 1
# sample estimates:
# HL2 of log(x^2) and log(y^2) 
#                   -0.1040422 
```

# Contributions

We are grateful for any contribution to the further development of the R package. If you experience any problems using the package or have suggestions for new features, please open an issue in the [issue tracker](https://github.com/s-abbas/robTests/issues). 

Authors
-------

**Sermad Abbas** ( [s-abbas](https://github.com/s-abbas) ) - *TU Dortmund University, Dortmund, Germany*, and 
**Barbara Brune** ( [b-brune](https://github.com/b-brune) ) - *TU Wien, Vienna, Austria*
