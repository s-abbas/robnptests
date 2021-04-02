
robTests
========

The R-package robTests contains different functions for robust two-sample tests for location.

Installation
------------

To install the package, the devtools package is required.

``` r
install.packages("devtools")
devtools::install_github("s-abbas/robTests")

library(robTests)
```

Use of the package
------------------

The robTests package contains various robust estimators for location and scale of a numeric sample, as well as a number of robust tests for the two sample location problem.

The location estimators implemented are:

-   the trimmed mean <code>trim\_mean()</code>
-   different asymmetrically trimmed means <code>asym\_trimmed\_mean()</code>
-   the winsorized mean <code>win\_mean()</code>
-   the one-sample Hodges-Lehmann estimator <code>hodges\_lehmann()</code> and the two-sample Hodges-Lehmann estimator for shift <code>hodges\_lehmann\_2sample()</code>
-   M-estimators with different tuning functions <code>m\_est()</code>

According scale estimators are:

-   the winsorized variance <code>win\_var()</code> and asymmetrically winsorized versions <code>asym\_win\_var()</code>
-   robust estimates for the joint scale of two samples based on median absolute deviations in <code>rob\_var()</code>

Based on these estimators, different robust and nonparametric tests for the two-sample location problem can be performed:

-   tests based on the median <code>med\_test()</code> and one- and two-sample Hodges-Lehmann estimators <code>hl1\_test()</code> and <code>hl2\_test()</code>
-   tests based on trimmed and asymmetrically trimmed means and variances <code>trimmed\_test()</code> and <code>asym\_trimmed\_test()</code>
-   hybrid tests using p-values derived from t- and trimmed t-statistics and test statistics based on Huber's M-estimator: <code>min\_c\_test()</code>, <code>min\_tc\_test()</code> and <code>min\_t\_test()</code>
-   permutation tests based on different M-estimators <code>m\_estimator\_tests()</code>

For details and references see the documentation of the different functions.

### Perform two sample tests for location shift

``` r
x <- rnorm(50); y <- rnorm(50)

## Perform a trimmed two sample test:
trimmed_test(x, y, gamma = 0.1)
#> 
#>  Trimmed two-sample t-test
#> 
#> data:  x and y
#> trimmed t = -1.7672, df = 54, p-value = 0.08285
#> alternative hypothesis: true difference in means is not equal to 0
#> sample estimates:
#> Trimmed mean of x Trimmed mean of y 
#>        -0.1805498         0.3931032

## Perform a radomization test based on the two sample Hodges-Lehmann estimator
hl2_test(x, y, method = "sampled", n.rep = 1000)
#> 
#>  Randomization test based on the Two-Sample Hodges-Lehmann
#>  estimator
#> 
#> data:  x and y
#> D = -0.53863, p-value = 0.06494
#> alternative hypothesis: true location shift is not equal to 0
#> sample estimates:
#> HL2 of x and y 
#>     -0.5977219
```

### Compute robust estimates of scale and location

``` r
hodges_lehmann(x)
#> [1] -0.1770372

trim_mean(x)
#> [1] -0.2028446

win_var(x)$var
#> [1] 1.099881
```

# Contributions

We are grateful for any contribution to the further development of the R package. If you experience any problems using the package or have suggestions for new features, please open an issue in the [issue tracker](https://github.com/s-abbas/robTests/issues). 

Authors
-------

**Sermad Abbas** ( [s-abbas](https://github.com/s-abbas) ) and **Barbara Brune** ( [b-brune](https://github.com/b-brune) ) - *TU Dortmund University, Germany*
