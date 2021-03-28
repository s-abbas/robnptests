
## Trimmed mean ----
testthat::test_that("trimmed_mean works correctly", {

  ## Generate exemplary input vector
  x <- c(1, 3, 8, 9, 15)

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

  # hodges_lehmann_2sample(x, y) should be equal to -hodges_lehmann_2sample(y, x)
  testthat::expect_equal(hodges_lehmann_2sample(x, y), -hodges_lehmann_2sample(y, x))

  # The two-sample Wilcoxon rank-sum statistic of aligned sampled needs to be
  # equal to its expected value under H_0 of two-sample Wilcoxon rank-sum test

  # Align second sample
  x <- x - hodges_lehmann_2sample(x, y)

  testthat::expect_equal(sum(rank(c(x, y))[1:m]), m * (m + n + 1)/2)
})

## M-estimator ----
testthat::test_that("m_est works correctly", {

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

