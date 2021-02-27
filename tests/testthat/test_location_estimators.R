
## Trimmed mean ----
testthat::test_that("trimmed_mean works correctly", {

  ## Generate exemplary input vector
  x <- c(1, 3, 8, 9, 15)

  ## The input arguments need to be checked

  # Missing input arguments
  testthat::expect_error(trim_mean(gamma = 0.2, na.rm = FALSE))
  testthat::expect_equal(trim_mean(x = x, na.rm = FALSE), trim_mean(x = x, gamma = 0.2, na.rm = FALSE))
  testthat::expect_equal(trim_mean(x = x, gamma = 0.2), trim_mean(x = x, gamma = 0.2, na.rm = FALSE))

  # Wrong data types
  testthat::expect_error(trim_mean(x = as.character(x), gamma = 0.2, na.rm = FALSE))
  testthat::expect_error(trim_mean(x = x, gamma = as.character(0.2), na.rm = FALSE))
  testthat::expect_error(trim_mean(x = x, gamma = 0.2, na.rm = 1))
  testthat::expect_error(trim_mean(x = x, gamma = 0.7, na.rm = FALSE), regexp = "'gamma' has to be in [0, 0.5].", fixed = TRUE)
  testthat::expect_error(trim_mean(x = x, gamma = -0.01, na.rm = FALSE), regexp = "'gamma' has to be in [0, 0.5].", fixed = TRUE)

  ## Missing values should be removed correctly
  testthat::expect_equal(trim_mean(x = c(x, NA), gamma = 0.2, na.rm = TRUE),
                         trim_mean(x = x, gamma = 0.2, na.rm = TRUE)
  )
  testthat::expect_equal(trim_mean(x = c(x, NA), gamma = 0.2, na.rm = FALSE), NA_real_)
  testthat::expect_equal(trim_mean(x = c(x, NA), gamma = 0.2),
                         trim_mean(x = c(x, NA), gamma = 0.2, na.rm = FALSE)
  )

  ## The computed value need to be identical to the output of mean(x, trim = ...)
  testthat::expect_equal(trim_mean(x = x, gamma = 0.2),
                         mean(x = x, trim = 0.2)
  )
  testthat::expect_equal(trim_mean(x = x, gamma = 0),
                         mean(x = x, trim = 0)
  )
  testthat::expect_equal(trim_mean(x = x, gamma = 0.5),
                         mean(x = x, trim = 0.5)
  )
})



## One-sample Hodges-Lehmann estimator
testthat::test_that("hodges_lehmann works correctly", {
  ## ___________________________________________________________________________
  ## Computation of one-sample Hodges-Lehmann estimator
  ## ___________________________________________________________________________
  x <- c(1, 3, 8, 9, 15)
  testthat::expect_equal(hodges_lehmann(x), 7)

  ## ___________________________________________________________________________
  ## Handling of missing values
  ## ___________________________________________________________________________
  y <- c(x, NA)
  testthat::expect_equal(hodges_lehmann(y, na.rm = TRUE), 7)
  testthat::expect_true(is.na(hodges_lehmann(y, na.rm = FALSE)))

  ## ___________________________________________________________________________
  ## Non-numeric input
  ## ___________________________________________________________________________
  z <- letters
  testthat::expect_error(hodges_lehmann(z))
}
)

## Two-sample Hodges-Lehmann estimator
testthat::test_that("hodges_lehmann_2sample works correctly", {

  testthat::skip_on_cran()

  ## ___________________________________________________________________________
  ## Computation of two-sample Hodges-Lehmann estimator
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(10)
  m <- length(x)
  y <- rnorm(12)
  n <- length(y)

  ## Align second sample
  x <- x - hodges_lehmann_2sample(x, y)

  ## Two-sample Wilcoxon rank-sum statistic of aligned sampled needs to be
  ## equal to its expected value under H_0 of two-sample Wilcoxon rank-sum test
  testthat::expect_equal(sum(rank(c(x, y))[1:m]), m * (m + n + 1)/2)

  ## ___________________________________________________________________________
  ## Handling of missing values
  ## ___________________________________________________________________________

  set.seed(108)

  x <- rnorm(10)
  m <- length(x)
  y <- c(rnorm(11), NA)
  n <- length(y)

  ## Align second sample after removing missing values
  x <- x - hodges_lehmann_2sample(x, y, na.rm = TRUE)

  ## Two-sample Wilcoxon rank-sum statistic of aligned samples needs to be
  ## equal to its expected value under H_0 of two-sample Wilcoxon rank-sum test
  testthat::expect_equal(sum(rank(c(x, y))[1:m]), m * (m + (n - 1) + 1)/2)

  ## If missing values are not removed, the function should return NA
  testthat::expect_true(is.na(hodges_lehmann_2sample(x, y, na.rm = FALSE)))

  ## ___________________________________________________________________________
  ## Non-numeric input
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(10)
  y <- letters[1:12]

  ## If there are non-numeric input values, the function should throw an error
  testthat::expect_error(hodges_lehmann_2sample(x, y))
}
)
