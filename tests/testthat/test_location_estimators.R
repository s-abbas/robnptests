
## Trimmed mean ----
testthat::test_that("trimmed_mean works correctly", {

  ## Generate exemplary input vector
  x <- c(1, 3, 8, 9, 15)

  ## The input arguments need to be checked

  # Missing input arguments
  testthat::expect_error(trim_mean(gamma = 0.2, na.rm = FALSE), regexp = "'x' is missing.", fixed = TRUE)
  testthat::expect_equal(trim_mean(x = x, na.rm = FALSE), trim_mean(x = x, gamma = 0.2, na.rm = FALSE))
  testthat::expect_equal(trim_mean(x = x, gamma = 0.2), trim_mean(x = x, gamma = 0.2, na.rm = FALSE))

  # Input arguments are NULL
  testthat::expect_error(trim_mean(x = NULL, gamma = 0.2, na.rm = FALSE), regexp = "'x' must not be NULL.", fixed = TRUE)
  testthat::expect_error(trim_mean(x = x, gamma = NULL, na.rm = FALSE), regexp = "'gamma' must not be NULL.", fixed = TRUE)
  testthat::expect_error(trim_mean(x = x, gamma = 0.2, na.rm = NULL), regexp = "'na.rm' must not be NULL.", fixed = TRUE)

  # Wrong data types
  testthat::expect_error(trim_mean(x = as.character(x), gamma = 0.2, na.rm = FALSE), regexp = "'x' has to a numeric vector.", fixed = TRUE)
  testthat::expect_error(trim_mean(x = x, gamma = as.character(0.2), na.rm = FALSE), regexp = "'gamma' has to be a numeric value.", fixed = TRUE)
  testthat::expect_error(trim_mean(x = x, gamma = 0.2, na.rm = 1), regexp = "'na.rm' has to be a logical value.", fixed = TRUE)

  # Length of input arguments
  testthat::expect_error(trim_mean(x = x, gamma = rep(0.1, 3)), regexp = "'gamma' has to be a single value, not a vector of length >= 1.", fixed = TRUE)
  testthat::expect_error(trim_mean(x = x, gamma = 0.2, na.rm = rep(TRUE, 3)), regexp = "'na.rm' has to be a single value, not a vector of length >= 1.", fixed = TRUE)

  # Wrong value for 'gamma'
  testthat::expect_error(trim_mean(x = x, gamma = 0.7, na.rm = FALSE), regexp = "'gamma' has to be a numeric value in [0, 0.5].", fixed = TRUE)
  testthat::expect_error(trim_mean(x = x, gamma = -0.01, na.rm = FALSE), regexp = "'gamma' has to be a numeric value in [0, 0.5].", fixed = TRUE)

  ## Missing values should be removed correctly
  testthat::expect_equal(trim_mean(x = c(x, NA), gamma = 0.2, na.rm = TRUE),
                         trim_mean(x = x, gamma = 0.2, na.rm = TRUE)
  )
  testthat::expect_equal(trim_mean(x = c(x, NA), gamma = 0.2, na.rm = FALSE), NA_real_)
  testthat::expect_equal(trim_mean(x = c(x, NA), gamma = 0.2),
                         trim_mean(x = c(x, NA), gamma = 0.2, na.rm = FALSE)
  )

  ## The computed value has to be identical to the output of mean(x, trim = ...)
  testthat::expect_identical(trim_mean(x = x, gamma = 0.2),
                             mean(x = x, trim = 0.2)
  )
  testthat::expect_identical(trim_mean(x = x, gamma = 0),
                             mean(x = x, trim = 0)
  )
  testthat::expect_identical(trim_mean(x = x, gamma = 0.5),
                             mean(x = x, trim = 0.5)
  )

  ## If all elements of 'x' are NA, the output needs to be NA
  testthat::expect_identical(trim_mean(x = NA_real_, gamma = 0.2, na.rm = TRUE), NA_real_)
  testthat::expect_identical(trim_mean(x = NA_real_, gamma = 0.2, na.rm = FALSE), NA_real_)
})

## Winsorized mean ----
testthat::test_that("win_mean works correctly", {

  ## Generate exemplary input vector
  x <- c(92, 19, 101, 58, 1053, 91, 26, 78, 10, 13, -40, 101, 86, 85, 15, 89, 89, 28, -5, 41)

  ## The input arguments need to be checked

  # Missing input arguments
  testthat::expect_error(win_mean(gamma = 0.2, na.rm = FALSE), regexp = "'x' is missing.", fixed = TRUE)
  testthat::expect_equal(win_mean(x = x, na.rm = FALSE), win_mean(x = x, gamma = 0.2, na.rm = FALSE))
  testthat::expect_equal(win_mean(x = x, gamma = 0.2), win_mean(x = x, gamma = 0.2, na.rm = FALSE))

  # Input arguments are NULL
  testthat::expect_error(win_mean(x = NULL, gamma = 0.2, na.rm = FALSE), regexp = "'x' must not be NULL.", fixed = TRUE)
  testthat::expect_error(win_mean(x = x, gamma = NULL, na.rm = FALSE), regexp = "'gamma' must not be NULL.", fixed = TRUE)
  testthat::expect_error(win_mean(x = x, gamma = 0.2, na.rm = NULL), regexp = "'na.rm' must not be NULL.", fixed = TRUE)

  # Wrong data types
  testthat::expect_error(win_mean(x = as.character(x), gamma = 0.2, na.rm = FALSE), regexp = "'x' has to a numeric vector.", fixed = TRUE)
  testthat::expect_error(win_mean(x = x, gamma = as.character(0.2), na.rm = FALSE), regexp = "'gamma' has to be a numeric value.", fixed = TRUE)
  testthat::expect_error(win_mean(x = x, gamma = 0.2, na.rm = 1), regexp = "'na.rm' has to be a logical value.", fixed = TRUE)

  # Length of input arguments
  testthat::expect_error(win_mean(x = x, gamma = rep(0.1, 3)), regexp = "'gamma' has to be a single value, not a vector of length >= 1.", fixed = TRUE)
  testthat::expect_error(win_mean(x = x, gamma = 0.2, na.rm = rep(TRUE, 3)), regexp = "'na.rm' has to be a single value, not a vector of length >= 1.", fixed = TRUE)

  # Wrong value for 'gamma'
  testthat::expect_error(win_mean(x = x, gamma = 0.7, na.rm = FALSE), regexp = "'gamma' has to be a numeric value in [0, 0.5].", fixed = TRUE)
  testthat::expect_error(win_mean(x = x, gamma = -0.01, na.rm = FALSE), regexp = "'gamma' has to be a numeric value in [0, 0.5].", fixed = TRUE)

  ## Missing values should be removed correctly
  testthat::expect_equal(win_mean(x = c(x, NA), gamma = 0.2, na.rm = TRUE),
                         win_mean(x = x, gamma = 0.2, na.rm = TRUE)
  )
  testthat::expect_equal(win_mean(x = c(x, NA), gamma = 0.2, na.rm = FALSE), NA_real_)
  testthat::expect_equal(win_mean(x = c(x, NA), gamma = 0.2),
                         win_mean(x = c(x, NA), gamma = 0.2, na.rm = FALSE)
  )

  ## The computed value has to be identical to 55.65 if 'gamma == 0.05'
  ## and to the sample mean if 'gamma == 0'
  testthat::expect_equal(win_mean(x = x, gamma = 0.05), 55.65)
  testthat::expect_equal(win_mean(x = x, gamma = 0), mean(x = x))

  ## If all elements of 'x' are NA, the output needs to be NA
  testthat::expect_identical(win_mean(x = NA_real_, gamma = 0.2, na.rm = TRUE), NA_real_)
  testthat::expect_identical(win_mean(x = NA_real_, gamma = 0.2, na.rm = FALSE), NA_real_)
})

## One-sample Hodges-Lehmann estimator ----
testthat::test_that("hodges_lehmann works correctly", {

  ## Generate exemplary input vector
  x <- c(1, 3, 8, 9, 15)

  ## The input arguments need to be checked

  # Missing input arguments
  testthat::expect_error(hodges_lehmann(na.rm = FALSE), regexp = "'x' is missing.", fixed = TRUE)

  # Input arguments are NULL
  testthat::expect_error(hodges_lehmann(x = NULL, na.rm = FALSE), regexp = "'x' must not be NULL.", fixed = TRUE)
  testthat::expect_error(hodges_lehmann(x = x, na.rm = NULL), regexp = "'na.rm' must not be NULL.", fixed = TRUE)

  # Wrong data types
  testthat::expect_error(hodges_lehmann(x = as.character(x), na.rm = FALSE), regexp = "'x' has to a numeric vector.", fixed = TRUE)
  testthat::expect_error(hodges_lehmann(x = x, na.rm = 1), regexp = "'na.rm' has to be a logical value.", fixed = TRUE)

  # Length of input arguments
  testthat::expect_error(hodges_lehmann(x = x[1], na.rm = FALSE), regexp = "'x' needs at least 2 values.", fixed = TRUE)
  testthat::expect_error(hodges_lehmann(x = x, na.rm = rep(TRUE, 3)), regexp = "'na.rm' has to be a single value, not a vector of length >= 1.", fixed = TRUE)

  ## Missing values should be removed correctly
  testthat::expect_equal(hodges_lehmann(x = c(x, NA), na.rm = TRUE),
                         hodges_lehmann(x = x, na.rm = TRUE)
  )
  testthat::expect_equal(hodges_lehmann(x = c(x, NA), na.rm = FALSE), NA_real_)
  testthat::expect_equal(hodges_lehmann(x = c(x, NA)), hodges_lehmann(x = c(x, NA), na.rm = FALSE)
  )

  ## The computed value has to be identical to 7
  testthat::expect_equal(hodges_lehmann(x), 7)

  ## If all elements of 'x' are NA, the output needs to be NA
  testthat::expect_identical(hodges_lehmann(x = rep(NA_real_, 2, na.rm = TRUE)), NA_real_)
  testthat::expect_identical(hodges_lehmann(x = rep(NA_real_, 2, na.rm = FALSE)), NA_real_)
})

## Two-sample Hodges-Lehmann estimator ----
testthat::test_that("hodges_lehmann_2sample works correctly", {

  ## Generate exemplary input vectors
  x <- c(1, 2, 5, 10, 14, 16, 7, 3, 9, 21)
  m <- length(x)
  y <- c(32, 3, 4, 8, 13, 12, 17, 20, 3, 11)
  n <- length(y)

  ## The input arguments need to be checked

  # Missing input arguments
  testthat::expect_error(hodges_lehmann_2sample(y = y, na.rm = FALSE), regexp = "'x' is missing.", fixed = TRUE)
  testthat::expect_error(hodges_lehmann_2sample(x = x, na.rm = FALSE), regexp = "'y' is missing.", fixed = TRUE)

  # Input arguments are NULL
  testthat::expect_error(hodges_lehmann_2sample(x = NULL, y = y, na.rm = FALSE), regexp = "'x' must not be NULL.", fixed = TRUE)
  testthat::expect_error(hodges_lehmann_2sample(x = x, y = NULL, na.rm = FALSE), regexp = "'y' must not be NULL.", fixed = TRUE)
  testthat::expect_error(hodges_lehmann_2sample(x = x, y = y, na.rm = NULL), regexp = "'na.rm' must not be NULL.", fixed = TRUE)

  # Wrong data types
  testthat::expect_error(hodges_lehmann_2sample(x = as.character(x), y = y, na.rm = FALSE), regexp = "'x' has to a numeric vector.", fixed = TRUE)
  testthat::expect_error(hodges_lehmann_2sample(x = x, y = as.character(y), na.rm = FALSE), regexp = "'y' has to a numeric vector.", fixed = TRUE)
  testthat::expect_error(hodges_lehmann_2sample(x = x, y = y, na.rm = 1), regexp = "'na.rm' has to be a logical value.", fixed = TRUE)

  # Length of input arguments
  testthat::expect_error(hodges_lehmann_2sample(x = x, y = y, na.rm = rep(TRUE, 3)), regexp = "'na.rm' has to be a single value, not a vector of length >= 1.", fixed = TRUE)

  ## Missing values should be removed correctly
  testthat::expect_equal(hodges_lehmann_2sample(x = c(x, NA), y = y, na.rm = TRUE),
                         hodges_lehmann_2sample(x = x, y = y, na.rm = TRUE)
  )
  testthat::expect_equal(hodges_lehmann_2sample(x = x, y = c(y, NA), na.rm = TRUE),
                         hodges_lehmann_2sample(x = x, y = y, na.rm = TRUE)
  )
  testthat::expect_equal(hodges_lehmann_2sample(x = c(x, NA), y = y, na.rm = FALSE), NA_real_)
  testthat::expect_equal(hodges_lehmann_2sample(x = x, y = c(y, NA), na.rm = FALSE), NA_real_)
  testthat::expect_equal(hodges_lehmann_2sample(x = c(x, NA), y = y), hodges_lehmann_2sample(x = c(x, NA), y = y, na.rm = FALSE))
  testthat::expect_equal(hodges_lehmann_2sample(x = x, y = c(y, NA)), hodges_lehmann_2sample(x = c(x, NA), y = y, na.rm = FALSE))

  ## Computation of two-sample Hodges-Lehmann estimator

  # Align second sample
  x <- x - hodges_lehmann_2sample(x, y)

  # The two-sample Wilcoxon rank-sum statistic of aligned sampled needs to be
  # equal to its expected value under H_0 of two-sample Wilcoxon rank-sum test
  testthat::expect_equal(sum(rank(c(x, y))[1:m]), m * (m + n + 1)/2)
})

## M-estimator ----
testthat::test_that("m_est works correctly", {

  ## The input arguments need to be checked

  # Generate exemplary input vector
  set.seed(108)
  x <- rnorm(10)

  # Missing input arguments
  testthat::expect_error(m_est(psi = "huber"), regexp = "'x' is missing.", fixed = TRUE)
  testthat::expect_error(m_est(x = x), regexp = "'psi' is missing.", fixed = TRUE)

  # Input arguments are NULL
  testthat::expect_error(m_est(x = NULL, psi = "huber"), regexp = "'x' must not be NULL", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = NULL), regexp = "'psi' must not be NULL", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", k = NULL), regexp = "'k' must not be NULL", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", tol = NULL), regexp = "'tol' must not be NULL", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", max.it = NULL), regexp = "'max.it' must not be NULL", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", na.rm = NULL), regexp = "'na.rm' must not be NULL", fixed = TRUE)

  # Input arguments are NA
  testthat::expect_error(m_est(x = x, psi = NA)) # An error is thrown by .Mpsi.tuning.default(psi)
  testthat::expect_error(m_est(x = x, psi = "huber", k = NA), regexp = "'k' must not be NA", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", tol = NA), regexp = "'tol' must not be NA", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", max.it = NA), regexp = "'max.it' must not be NA", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", na.rm = NA), regexp = "'na.rm' must not be NA", fixed = TRUE)

  # Wrong data types
  testthat::expect_error(m_est(x = as.character(x), psi = "huber"), regexp = "'x' has to a numeric vector.", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = 1), regexp = "'psi' has to a character value.", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", k = "1.345"), regexp = "'k' has to be a numeric value.", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", tol = "1e-06"), regexp = "'tol' has to be a numeric value.", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", max.it = "15"), regexp = "'max.it' has to be a numeric value.", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", na.rm = 1), regexp = "'na.rm' has to be a logical value.", fixed = TRUE)

  # Length of input arguments
  testthat::expect_error(m_est(x = x, psi = c("huber", "hampel", "bisquare"))) # An error is thrown by .Mpsi.tuning.defaul(psi)
  testthat::expect_error(m_est(x = x, psi = "welsh"), regexp = "'psi' must be one of 'huber', 'hampel', or 'bisquare'.", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = c("tukey"))) # An error is thrown by .Mpsi.tuning.defaul(psi)
  testthat::expect_error(m_est(x = x, psi = "huber", k = c(1.03, 1.345)), regexp = "'k' has to be a single value, not a vector of length >= 1.", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", tol = c(1e-06, 1e-10)), regexp = "'tol' has to be a single value, not a vector of length >= 1.", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", max.it = c(10, 15)), regexp = "'max.it' has to be a single value, not a vector of length >= 1.", fixed = TRUE)
  testthat::expect_error(m_est(x = x, psi = "huber", na.rm = c(TRUE, FALSE)), regexp = "'na.rm' has to be a single value, not a vector of length >= 1.", fixed = TRUE)

  ## Missing values should be removed correctly
  testthat::expect_equal(m_est(x = c(x, NA), psi = "huber", na.rm = TRUE),
                         m_est(x = x, psi = "huber")
  )
  testthat::expect_equal(m_est(x = c(x, NA), psi = "huber", na.rm = FALSE),
                         NA_real_
  )
  testthat::expect_equal(m_est(x = c(x, NA), psi = "huber"),
                         m_est(x = c(x, NA), psi = "huber", na.rm = FALSE)
  )

  ## TODO: Output needs to be tested.
})

