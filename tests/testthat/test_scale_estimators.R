# Winsorized variance ----
testthat::test_that("win_var works correctly", {

  # Exemplary input vectors ----

  # Even sample size
  x <- c(19, 72, 51, 43, 87, 91, 50, 38, 4, 63)

  # Odd sample size
  x1 <- c(19, 72, 51, 43, 87, 91, 50, 38, 4, 63, 75)

  # Removing missing values ----
  testthat::expect_equal(win_var(x = c(x, NA), na.rm = TRUE), win_var(x = x))

  testthat::expect_equal(win_var(x = c(x, NA)),
                         list(var = NA_real_, h = NA_real_))

  testthat::expect_equal(win_var(x = c(x, NA)),
                         win_var(x = c(x, NA), na.rm = FALSE))

  # All values in 'x' are equal ----

  # No NAs in 'x'
  testthat::expect_error(win_var(x = rep(0, 5)))

  # NAs in 'x'
  testthat::expect_error(win_var(x = c(NA, rep(0, 5)), na.rm = TRUE))

  # Location invariance, scale equivariance, and permutation invariance ----

  # Location invariance
  testthat::expect_equal(win_var(x = x + 5, gamma = 0),
                         win_var(x = x, gamma = 0))

  testthat::expect_equal(win_var(x = x + 5, gamma = 0.25),
                         win_var(x = x, gamma = 0.25))

  # Scale equivariance
  testthat::expect_equal(win_var(x = 2 * x, gamma = 0),
                         list(var = 4 * win_var(x, gamma = 0)$var,
                              h = length(x)))

  testthat::expect_equal(win_var(x = 2 * x, gamma = 0.25),
                         list(var = 4 * win_var(x, gamma = 0.25)$var, h = 6))

  # Permutation invariance
  testthat::expect_equal(win_var(x = x, gamma = 0),
                         win_var(x = sort(x), gamma = 0))

  testthat::expect_equal(win_var(x = x, gamma = 0.25),
                         win_var(x = sort(x), gamma = 0.25))

  # Check output ----

  # For 'gamma' = 0, the winsorized variance of 'x' is equal to the sample variance
  testthat::expect_equal(win_var(x = x), list(var = var(x), h = length(x)))

  # For 'gamma' = 0.2, the winsorized variance of 'x' is equal to 218.4556
  testthat::expect_equal(win_var(x = x, gamma = 0.2),
                         list(var = 218.4556, h = 6), tolerance = 1e-6)

  # For 'gamma' = 0.3, the winsorized variance of 'x' is equal to 90.05556
  testthat::expect_equal(win_var(x = x, gamma = 0.3),
                         list(var = 90.05556, h = 4), tolerance = 1e-6)

  # For 'gamma' = 0.25, the winsorized variance of 'x1' is equal to 90.05556
  testthat::expect_equal(win_var(x = x1, gamma = 0.25),
                         list(var = 258.9636, h = 7), tolerance = 1e-6)

  # The output should be a list of numeric scalars
  checkmate::expect_list(win_var(x = x), types = rep("numeric", 2))

  checkmate::expect_list(win_var(x = c(x, NA_real_), na.rm = FALSE),
                         types = rep("numeric", 2))
})

# Robust variance for 'hl1_test', 'hl2_test', and 'med_test' ----
testthat::test_that("rob_scale works correctly", {

  # Exemplary input vectors ----
  x <- c(7, 2, 1, 6, 8)
  y <- c(5, 9, 6, 7, 8)

  # Removing missing values ----
  testthat::expect_equal(rob_scale(x = c(x, NA), y = c(y, NA), na.rm = TRUE),
                         rob_scale(x = x, y = y))

  testthat::expect_equal(rob_scale(x = c(x, NA), y = c(y, NA)), NA_real_)

  testthat::expect_equal(rob_scale(x = c(x, NA), y = c(y, NA)),
                         rob_scale(x = c(x, NA), y = c(y, NA), na.rm = FALSE))

  # All values in 'x' and 'y' are equal ----

  # No NAs in 'x' and 'y'
  testthat::expect_equal(rob_scale(x = rep(0, 5), y = rep(0, 5)), 0)
  testthat::expect_error(rob_scale(x = rep(0, 5), y = rep(0, 5),
                                   check.for.zero = TRUE))

  # NAs in 'x' and 'y'
  testthat::expect_equal(rob_scale(x = c(NA, rep(0, 5)), y = c(NA, rep(0, 5)),
                                 na.rm = TRUE), 0)

  testthat::expect_error(rob_scale(x = c(NA, rep(0, 5)), y = c(NA, rep(0, 5)),
                                 na.rm = TRUE, check.for.zero = TRUE))

  # Location invariance, scale equivariance, and permutation invariance ----
  types <- c("S1", "S2", "S3", "S4")

  for (i in seq_along(types)) {
    # Location invariance
    testthat::expect_equal(rob_scale(x = x + 2, y = y + 3, type = types[i]),
                           rob_scale(x = x, y = y, type = types[i]))

    # Scale equivariance
    testthat::expect_equal(rob_scale(x = 2 * x, y = 2 * y, type = types[i]),
                           2 * rob_scale(x = x, y = y, type = types[i]))

    # Permutation invariance
    testthat::expect_equal(rob_scale(x = sort(x), y = sort(y), type = types[i]),
                           rob_scale(x = x, y = y, type = types[i]))

    # Switch 'x' and 'y'
    testthat::expect_equal(rob_scale(x = y, y = x, type = types[i]),
                           rob_scale(x = x, y = y, type = types[i]))
  }

  # Check output ----

  # For 'type' = 'S1', the output should be 2
  testthat::expect_equal(rob_scale(x = x, y = y, type = "S1"), 2)

  # For 'type' = 'S2', the output should be 2
  testthat::expect_equal(rob_scale(x = x, y = y, type = "S2"), 2)

  # For 'type' = 'S3', the output should be 3
  testthat::expect_equal(rob_scale(x = x, y = y, type = "S3"), 3)

  # For 'type' = 'S4', the output should be 3
  testthat::expect_equal(rob_scale(x = x, y = y, type = "S4"), 3)

  for (i in seq_along(types)) {
    # The output should be a numeric scalar
    checkmate::expect_numeric(rob_scale(x = x, y = y, type = types[i]), len = 1)

    checkmate::expect_numeric(rob_scale(x = c(x, NA_real_), y = y,
                                      type = types[i]),
                              len = 1)
  }

  # An error is expected when the scale estimate is zero even though the data
  # are not constant
  x1 <- c(7, 2, 1, 2, 2)
  y1 <- c(5, 9, 2, 2, 2)
  testthat::expect_error(rob_scale(x = x1, y = y1, type = "S3", check.for.zero = TRUE))
})
