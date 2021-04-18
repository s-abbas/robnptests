## 'preprocess_data' ----
testthat::test_that("preprocess_data works correctly", {

  ## Generate exemplary input vectors
  set.seed(108)
  x <- rnorm(10)
  y <- rnorm(10)

  ## Removing missing values ----

  # No missing values
  testthat::expect_equal(preprocess_data(x = x, y = y, delta = 0, na.rm = TRUE,
                                         wobble = FALSE, var.test = FALSE),
                         list(x = x, y = y, delta = 0)
  )

  # Missing values and remove missing values
  testthat::expect_equal(preprocess_data(x = c(x, NA_real_), y = c(y, NA_real_),
                                         delta = 0, na.rm = TRUE, wobble = FALSE,
                                         var.test = FALSE),
                         list(x = x, y = y, delta = 0)
  )

  # Missing values and do not remove missing values
  testthat::expect_equal(preprocess_data(x = c(x, NA_real_), y = c(y, NA_real_),
                                         delta = 0, na.rm = FALSE, wobble = FALSE,
                                         var.test = FALSE),
                         list(x = c(x, NA_real_), y = c(y, NA_real_), delta = 0)
  )

  ## Warning, if the samples contain less than five observations ----

  # No missing values
  testthat::expect_error(preprocess_data(x = x[1:4], y = y, delta = 0, na.rm = TRUE,
                                         wobble = FALSE, var.test = FALSE)
  )

  # Too few observations after removing missing values
  testthat::expect_error(preprocess_data(x = c(x[1:4], rep(NA_real_, 2)), y = y,
                                         delta = 0, na.rm = TRUE,
                                         wobble = FALSE, var.test = FALSE)
  )

  ## Wobbling ----

  # No duplicated values in input vectors -> no wobbling
  testthat::expect_equal(preprocess_data(x = c(x, NA_real_), y = c(y, NA_real_),
                                         delta = 0, na.rm = TRUE, wobble = TRUE,
                                         var.test = FALSE),
                         list(x = x, y = y, delta = 0)
  )

  # Duplicated values within each sample
  testthat::expect_warning(preprocess_data(x = c(x, x[10]), y = c(y, y[10]),
                                           delta = 0, na.rm = TRUE, wobble = TRUE,
                                           wobble.seed = 123, var.test = FALSE)
  )

  # Duplicated values between the samples
  testthat::expect_warning(preprocess_data(x = c(x, x[10]), y = c(y, x[10]),
                                           delta = 0, na.rm = TRUE, wobble = TRUE,
                                           wobble.seed = 123, var.test = FALSE)
  )

  ## Transformation to test for difference in scale ----

  # No zeros in samples
  testthat::expect_equal(preprocess_data(x = x, y = y, delta = 1, na.rm = TRUE,
                                         wobble = FALSE, var.test = TRUE),
                         list(x = log(x^2), y = log(y^2), delta = 0)
  )

  # No zeros in samples but 'delta' = 0
  testthat::expect_error(preprocess_data(x = x, y = y, delta = 0, na.rm = TRUE,
                                         wobble = FALSE, var.test = TRUE)
  )

  # Zeros in samples and 'delta' != 0
  testthat::expect_warning(preprocess_data(x = c(x, 0), y = c(y, 0), delta = 1,
                                           na.rm = TRUE, wobble = FALSE,
                                           wobble.seed = 123, var.test = TRUE)
  )
})

## 'select_method' ----
testthat::test_that("select_method works correctly", {

  ## Generate exemplary input vectors

  # Sample sizes m = n = 10
  set.seed(108)
  x1 <- rnorm(10)
  y1 <- rnorm(10)

  # Sample sizes m = n = 30
  set.seed(108)
  x2 <- rnorm(30)
  y2 <- rnorm(30)

  ## Principle specified by user ----
  testthat::expect_equal(select_method(x = x1, y = y1, method = "randomization",
                                       test.name = "hl1_test", n.rep = 10000),
                         "randomization")

  ## Automatic selection of the principle ----

  # Automatic selection is implemented for the test
  # m = n = 10 -> randomization test
  testthat::expect_equal(select_method(x = x1, y = y1,
                                       method = c("asymptotic", "permutation",
                                                  "randomization"),
                                       test.name = "hl1_test", n.rep = 10000),
                         "randomization")

  # m = n = 5 and n.rep >= choose(m+n, n) -> permutation test
  testthat::expect_equal(select_method(x = x1, y = y1,
                                      method = c("asymptotic", "permutation",
                                                 "randomization"),
                                      test.name = "hl1_test", n.rep = 200000),
                        "permutation")

  # m = n = 30 -> asymptotic test
  testthat::expect_equal(select_method(x = x2, y = y2,
                                       method = c("asymptotic", "permutation",
                                                  "randomization"),
                                       test.name = "hl1_test", n.rep = 10000),
                         "asymptotic")

  # Automatic selection is not implemented for the test
  testthat::expect_error(select_method(x = x1, y = y1,
                                       method = c("asymptotic", "permutation",
                                                  "randomization"),
                                       test.name = "m_test1", n.rep = 10000)
  )


})
