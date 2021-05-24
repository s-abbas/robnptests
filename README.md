
robnptests
========

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/s-abbas/robnptests/branch/master/graph/badge.svg)](https://codecov.io/gh/s-abbas/robnptests?branch=develop)
<!-- badges: end -->

The R package `robnptests` contains different robust and non-parametric tests for the two-sample location problem. The tests allow for comparisons of either the location or the scale parameters of two random samples.

Installation
------------

To install the package, the `devtools` package is required.

``` r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github("s-abbas/robnptests")

library(robnptests)
```

Scope and Usage
---------------

The robust and non-parametric tests in this R package follow the construction principle of the popular t-test: A robust estimate for the location difference of the two samples is divided by a robust estimate of scale. 
The p-values can either be computed using the permutation principle, the randomization principle, or the asymptotic distribution of the estimators. 
If the principle to compute the p-value is not specified by the user, it will be selected automatically depending on the sample size. 
The functions used to compute the location and scale estimates are also made available to the user.

The following list shows the currently implemented tests in the package:

* tests based on the median (`med_test`), the one-sample Hodges-Lehmann estimator (`hl1_test`), and the two-sample Hodges-Lehmann estimator (`hl2_test`), scaled by  robust estimators
* Yuen's t-test (`trimmed_test`) 
* tests based on the Huber-, Hampel- or Bisquare-M-estimator (`m_test`).

Even though the test statistics compare location estimates of the samples, they can be used to identify scale differences.
This is achieved by setting the argument `var.test = TRUE`, with which the observations in the samples are log-transformed so that scale differences between the original samples correspond to location differences in the transformed samples.

Details on the tests and references can be found on the help pages of the functions and the vignette `vignette(robnptests)`.

### Example 1: Asymptotic test for location difference using the two-sample Hodges-Lehmann estimator

``` r
set.seed(121)
x <- rnorm(50); y <- rnorm(50)

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
```

### Example 2: Asymptotic test for scale difference using the two-sample Hodges-Lehmann estimator

``` r
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

We are grateful for any contribution to the further development of the R package. If you experience any problems using the package or have suggestions for new features, please open an issue in the [issue tracker](https://github.com/s-abbas/robnptests/issues). 

Authors
-------

**Sermad Abbas** ( [s-abbas](https://github.com/s-abbas) ) - *TU Dortmund University, Dortmund, Germany*, and 
**Barbara Brune** ( [b-brune](https://github.com/b-brune) ) - *TU Wien, Vienna, Austria*
