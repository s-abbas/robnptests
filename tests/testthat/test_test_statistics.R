## Trimmed t ----
testthat::test_that("trimmed_t works correctly", {

  ## Generate exemplary input vectors ----
  x <- c(22, 23, 31, 71, 74)
  y <- c(4, 4, 8, 19, 26, 74)

  # Results for 'x' computed by hand with 'gamma = 0.2':
  # trimmed_mean(x) = 41.66667
  # win_var(x) = 627.2
  # h = 3
  #
  # Results for 'y' computed by hand with 'gamma = 0.2':
  # trimmed_mean(y) = 14.25
  # win_var(y) = 109.5
  # h = 4
  #
  # Test statistic: 1.451922

  ## Removal of missing values ----
  testthat::expect_equal(trimmed_t(x = c(x, NA), y = y, na.rm = TRUE),
                         trimmed_t(x = x, y = y))

  testthat::expect_equal(trimmed_t(x = c(x, NA), y = y, na.rm = FALSE),
                         list(statistic = NA_real_,
                              estimates = c(NA_real_, trim_mean(y)),
                              df = NA_real_))

  ## Location and scale invariance of the test statistic ----
  testthat::expect_equal(trimmed_t(x = 2 * x + 3, y = 2 * y + 3)$statistic,
                         trimmed_t(x = x, y = y)$statistic)

  ## Output check ----

  # The output should be a list of numeric scalars
  checkmate::expect_list(trimmed_t(x = x, y = y), types = rep("numeric", 3))

  checkmate::expect_list(trimmed_t(x = c(NA, x), y = c(NA, y)),
                         types = rep("numeric", 3))

  # The output should be equal to the manually computed values
  testthat::expect_equal(trimmed_t(x = x, y = y, gamma = 0.2),
                         list(statistic = 1.451922,
                              estimates = c(41.66667, 14.25),
                              df = 5),
                         tol = 1e-07)

})
