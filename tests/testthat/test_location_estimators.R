# Trimmed mean ----
testthat::test_that("trimmed_mean works correctly", {

  # Exemplary input vector ----
  x <- c(1, 3, 8, 9, 15)

  # Removing missing values ----
  testthat::expect_equal(trim_mean(x = c(x, NA), gamma = 0.2, na.rm = TRUE),
                         trim_mean(x = x, gamma = 0.2, na.rm = TRUE)
  )
  testthat::expect_equal(trim_mean(x = c(x, NA), gamma = 0.2, na.rm = FALSE), NA_real_)
  testthat::expect_equal(trim_mean(x = c(x, NA), gamma = 0.2),
                         trim_mean(x = c(x, NA), gamma = 0.2, na.rm = FALSE)
  )

  # Check output ----

  # The computed value has to be identical to the output of mean(x, trim = ...)
  testthat::expect_identical(trim_mean(x = x, gamma = 0.2),
                             mean(x = x, trim = 0.2)
  )
  testthat::expect_identical(trim_mean(x = x, gamma = 0),
                             mean(x = x, trim = 0)
  )
  testthat::expect_identical(trim_mean(x = x, gamma = 0.5),
                             mean(x = x, trim = 0.5)
  )

  # The estimator should be location and scale equivariant
  testthat::expect_equal(trim_mean(x = 3 * x + 2), 3 * trim_mean(x = x) + 2)

  # The output should be a numeric scalar
  checkmate::expect_numeric(trim_mean(x = x), len = 1)
  checkmate::expect_numeric(trim_mean(x = c(x, NA), na.rm = FALSE), len = 1)
})

# Winsorized mean ----
testthat::test_that("win_mean works correctly", {

  # Exemplary input vector ----
  x <- c(92, 19, 101, 58, 1053, 91, 26, 78, 10, 13, -40, 101, 86, 85, 15, 89,
         89, 28, -5, 41)

  # Removing missing values ----
  testthat::expect_equal(win_mean(x = c(x, NA), gamma = 0.2, na.rm = TRUE),
                         win_mean(x = x, gamma = 0.2, na.rm = TRUE)
  )
  testthat::expect_equal(win_mean(x = c(x, NA), gamma = 0.2, na.rm = FALSE), NA_real_)
  testthat::expect_equal(win_mean(x = c(x, NA), gamma = 0.2),
                         win_mean(x = c(x, NA), gamma = 0.2, na.rm = FALSE)
  )

  # Check output ----

  # The computed value has to be identical to 55.65 if 'gamma' = 0.05
  # and to the sample mean if 'gamma' = 0
  testthat::expect_equal(win_mean(x = x, gamma = 0.05), 55.65)
  testthat::expect_equal(win_mean(x = x, gamma = 0), mean(x = x))

  # The estimator should be location and scale equivariant
  testthat::expect_equal(win_mean(x = 3 * x + 2), 3 * win_mean(x = x) + 2)

  # The output should be a numeric scalar
  checkmate::expect_numeric(win_mean(x = x), len = 1)
  checkmate::expect_numeric(win_mean(x = c(x, NA), na.rm = FALSE), len = 1)
})

# One-sample Hodges-Lehmann estimator ----
testthat::test_that("hodges_lehmann works correctly", {

  # Generate exemplary input vector ----
  x <- c(1, 3, 8, 9, 15)

  # Removing missing values ----
  testthat::expect_equal(hodges_lehmann(x = c(x, NA), na.rm = TRUE),
                         hodges_lehmann(x = x, na.rm = TRUE)
  )
  testthat::expect_equal(hodges_lehmann(x = c(x, NA), na.rm = FALSE), NA_real_)
  testthat::expect_equal(hodges_lehmann(x = c(x, NA)), hodges_lehmann(x = c(x, NA), na.rm = FALSE)
  )

  # Check output ----

  # The computed value has to be identical to 7
  testthat::expect_equal(hodges_lehmann(x), 7)

  # The estimator should be location and scale equivariant
  testthat::expect_equal(hodges_lehmann(x = 3 * x + 2), 3 * hodges_lehmann(x = x) + 2)

  # The output should be a numeric scalar
  checkmate::expect_numeric(hodges_lehmann(x = x), len = 1)
  checkmate::expect_numeric(hodges_lehmann(x = c(x, NA), na.rm = FALSE), len = 1)
})

# Two-sample Hodges-Lehmann estimator ----
testthat::test_that("hodges_lehmann_2sample works correctly", {

  # Generate exemplary input vectors ----
  x <- c(1, 2, 5, 10, 14, 16, 7, 3, 9, 21)
  m <- length(x)
  y <- c(32, 3, 4, 8, 13, 12, 17, 20, 3, 11)
  n <- length(y)

  # Removing missing values ----
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

  # Check output ----

  # hodges_lehmann_2sample(x, y) should be equal to -hodges_lehmann_2sample(y, x)
  testthat::expect_equal(hodges_lehmann_2sample(x, y), -hodges_lehmann_2sample(y, x))

  # The two-sample Wilcoxon rank-sum statistic of aligned sampled needs to be
  # equal to its expected value under H_0

  # Align second sample
  x <- x - hodges_lehmann_2sample(x, y)

  testthat::expect_equal(sum(rank(c(x, y))[1:m]), m * (m + n + 1)/2)

  # Estimator should be location invariant and scale equivariant
  testthat::expect_equal(hodges_lehmann_2sample(x = x + 2, y = y + 2),
                         hodges_lehmann_2sample(x = x, y = y))
  testthat::expect_equal(hodges_lehmann_2sample(x = 3 * x, y = 3 * y),
                         3 * hodges_lehmann_2sample(x = x, y = y))

  # The output should be a numeric scalar
  checkmate::expect_numeric(hodges_lehmann_2sample(x = x, y = y), len = 1)
  checkmate::expect_numeric(hodges_lehmann_2sample(x = c(x, NA), y = y, na.rm = FALSE), len = 1)
})

# M-estimator ----
testthat::test_that("m_est works correctly", {
  psi.funs <- c("huber", "hampel", "bisquare")

  # Exemplary input vector ----
  set.seed(108)
  x <- rnorm(10)

  for (i in seq_along(psi.funs)) {

    # Removing missing values ----
    testthat::expect_equal(m_est(x = c(x, NA), psi = psi.funs[i], na.rm = TRUE),
                           m_est(x = x, psi = psi.funs[i])
    )
    testthat::expect_equal(m_est(x = c(x, NA), psi = psi.funs[i], na.rm = FALSE),
                           NA_real_
    )
    testthat::expect_equal(m_est(x = c(x, NA), psi = psi.funs[i]),
                           m_est(x = c(x, NA), psi = psi.funs[i], na.rm = FALSE)
    )

    # Check output ----

    # Location estimators should be location and scale equivariant
    testthat::expect_equal(m_est(x = 3 * x + 2, psi = psi.funs[i])$est,
                           3 * m_est(x = x, psi = psi.funs[i])$est + 2, tolerance = 10e-6)

    # Variance estimators should be location invariant and scale equivariant
    testthat::expect_equal(m_est(x = 3 * x + 2, psi = psi.funs[i])$var,
                           9 * m_est(x = x, psi = psi.funs[i])$var, tolerance = 10e-6)

    # The output should be a list contain two numeric values
    checkmate::expect_list(m_est(x = x, psi = psi.funs[i]), types = rep("numeric", 2))

    # Create and compare snapshots of test output
    # The output of the function cannot be computed manually
    testthat::expect_snapshot_output(m_est(x = x, psi = psi.funs[i]))
  }
})
