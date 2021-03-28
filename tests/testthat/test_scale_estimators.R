
## Winsorized variance ----
testthat::test_that("win_var works correctly", {

  ## Generate exemplary input vector
  x <- c(19, 72, 51, 43, 87, 91, 50, 38, 4, 63)
  x1 <- c(19, 72, 51, 43, 87, 91, 50, 38, 4, 63, 75)

  ## Missing values should be removed correctly ----
  testthat::expect_equal(win_var(x = c(x, NA), na.rm = TRUE), win_var(x = x))
  testthat::expect_equal(win_var(x = c(x, NA)), list(var = NA_real_, h = NA_real_))
  testthat::expect_equal(win_var(x = c(x, NA)), win_var(x = c(x, NA), na.rm = FALSE))

  ## The estimator should be location and permutation invariant ----
  # Location invariance
  testthat::expect_equal(win_var(x = x + 5, gamma = 0), win_var(x = x, gamma = 0))
  testthat::expect_equal(win_var(x = x + 5, gamma = 0.25), win_var(x = x, gamma = 0.25))

  # Scale equivariance
  testthat::expect_equal(win_var(x = 2 * x, gamma = 0), list(var = 4 * win_var(x, gamma = 0)$var, h = 10))
  testthat::expect_equal(win_var(x = 2 * x, gamma = 0.25), list(var = 4 * win_var(x, gamma = 0.25)$var, h = 6))

  # Permutation invariance
  testthat::expect_equal(win_var(x = x, gamma = 0), win_var(x = sort(x), gamma = 0))
  testthat::expect_equal(win_var(x = x, gamma = 0.25), win_var(x = sort(x), gamma = 0.25))

  ## Check output ----

  # For 'gamma' = 0, the winsorized variance of 'x' is equal to the sample variance
  testthat::expect_equal(win_var(x = x), list(var = var(x), h = length(x)))

  # For 'gamma' = 0.2, the winsorized variance of 'x' is equal to 218.4556
  testthat::expect_equal(win_var(x = x, gamma = 0.2), list(var = 218.4556, h = 6), tol = 1e-6)

  # For 'gamma' = 0.3, the winsorized variance of 'x' is equal to 90.05556
  testthat::expect_equal(win_var(x = x, gamma = 0.3), list(var = 90.05556, h = 4), tol = 1e-6)

  # For 'gamma' = 0.25, the winsorized variance of 'x1' is equal to 90.05556
  testthat::expect_equal(win_var(x = x1, gamma = 0.25), list(var = 258.9636, h = 7), tol = 1e-6)

  # The output should be a list of numeric scalars
  checkmate::expect_list(win_var(x = x), types = rep("numeric", 2))
})

testthat::test_that("Robust scale estimators work correctly", {

  # testthat::skip_on_cran()

  ## Generate exemplary input vector
  x <- c(7, 2, 1, 6, 8)
  y <- c(5, 9, 6, 7, 8)

  ## Missing values should be removed correctly ----
  testthat::expect_equal(rob_var(x = c(x, NA), y = y, na.rm = TRUE), rob_var(x = x, y = y))
  testthat::expect_equal(rob_var(x = c(x, NA), y = y), NA_real_)
  testthat::expect_equal(rob_var(x = c(x, NA), y = y), rob_var(x = c(x, NA), y = y, na.rm = FALSE))

  testthat::expect_equal(rob_var(x = x, y = c(y, NA), na.rm = TRUE), rob_var(x = x, y = y))
  testthat::expect_equal(rob_var(x = x, y = c(y, NA)), NA_real_)
  testthat::expect_equal(rob_var(x = x, y = c(y, NA)), rob_var(x = x, y = c(y, NA), na.rm = FALSE))

  testthat::expect_equal(rob_var(x = c(x, NA), y = c(y, NA), na.rm = TRUE), rob_var(x = x, y = y))
  testthat::expect_equal(rob_var(x = c(x, NA), y = c(y, NA)), NA_real_)
  testthat::expect_equal(rob_var(x = c(x, NA), y = c(y, NA)), rob_var(x = c(x, NA), y = c(y, NA)), na.rm = FALSE)

  ## The estimator should be location and permutation invariant ----
  # Location invariance
  testthat::expect_equal(rob_var(x = x + 5, y = y, type = "S1"), rob_var(x = x, y = y, type = "S1"))
  testthat::expect_equal(rob_var(x = x + 5, y = y, type = "S2"), rob_var(x = x, y = y, type = "S2"))
  testthat::expect_equal(rob_var(x = x + 5, y = y, type = "S3"), rob_var(x = x, y = y, type = "S3"))
  testthat::expect_equal(rob_var(x = x + 5, y = y, type = "S4"), rob_var(x = x, y = y, type = "S4"))

  testthat::expect_equal(rob_var(x = x, y = y + 5, type = "S1"), rob_var(x = x, y = y, type = "S1"))
  testthat::expect_equal(rob_var(x = x, y = y + 5, type = "S2"), rob_var(x = x, y = y, type = "S2"))
  testthat::expect_equal(rob_var(x = x, y = y + 5, type = "S3"), rob_var(x = x, y = y, type = "S3"))
  testthat::expect_equal(rob_var(x = x, y = y + 5, type = "S4"), rob_var(x = x, y = y, type = "S4"))

  testthat::expect_equal(rob_var(x = x + 2, y = y + 3, type = "S1"), rob_var(x = x, y = y, type = "S1"))
  testthat::expect_equal(rob_var(x = x + 2, y = y + 3, type = "S2"), rob_var(x = x, y = y, type = "S2"))
  testthat::expect_equal(rob_var(x = x + 2, y = y + 3, type = "S3"), rob_var(x = x, y = y, type = "S3"))
  testthat::expect_equal(rob_var(x = x + 2, y = y + 3, type = "S4"), rob_var(x = x, y = y, type = "S4"))

  # Permutation invariance
  testthat::expect_equal(rob_var(x = sort(x), y = y, type = "S1"), rob_var(x = x, y = y, type = "S1"))
  testthat::expect_equal(rob_var(x = sort(x), y = y, type = "S2"), rob_var(x = x, y = y, type = "S2"))
  testthat::expect_equal(rob_var(x = sort(x), y = y, type = "S3"), rob_var(x = x, y = y, type = "S3"))
  testthat::expect_equal(rob_var(x = sort(x), y = y, type = "S4"), rob_var(x = x, y = y, type = "S4"))

  testthat::expect_equal(rob_var(x = x, y = sort(y), type = "S1"), rob_var(x = x, y = y, type = "S1"))
  testthat::expect_equal(rob_var(x = x, y = sort(y), type = "S2"), rob_var(x = x, y = y, type = "S2"))
  testthat::expect_equal(rob_var(x = x, y = sort(y), type = "S3"), rob_var(x = x, y = y, type = "S3"))
  testthat::expect_equal(rob_var(x = x, y = sort(y), type = "S4"), rob_var(x = x, y = y, type = "S4"))

  testthat::expect_equal(rob_var(x = sort(x), y = sort(y), type = "S1"), rob_var(x = x, y = y, type = "S1"))
  testthat::expect_equal(rob_var(x = sort(x), y = sort(y), type = "S2"), rob_var(x = x, y = y, type = "S2"))
  testthat::expect_equal(rob_var(x = sort(x), y = sort(y), type = "S3"), rob_var(x = x, y = y, type = "S3"))
  testthat::expect_equal(rob_var(x = sort(x), y = sort(y), type = "S4"), rob_var(x = x, y = y, type = "S4"))

  # Scale equivariance
  testthat::expect_equal(rob_var(x = 2 * x, y = 2 * y, type = "S1"), 2 * rob_var(x = x, y = y, type = "S1"))
  testthat::expect_equal(rob_var(x = 2 * x, y = 2 * y, type = "S2"), 2 * rob_var(x = x, y = y, type = "S2"))
  testthat::expect_equal(rob_var(x = 2 * x, y = 2 * y, type = "S3"), 2 * rob_var(x = x, y = y, type = "S3"))
  testthat::expect_equal(rob_var(x = 2 * x, y = 2 * y, type = "S4"), 2 * rob_var(x = x, y = y, type = "S4"))

  # Switch 'x' and 'y'
  testthat::expect_equal(rob_var(x = y, y = x, type = "S1"), rob_var(x = x, y = y, type = "S1"))
  testthat::expect_equal(rob_var(x = y, y = x, type = "S2"), rob_var(x = x, y = y, type = "S2"))
  testthat::expect_equal(rob_var(x = y, y = x, type = "S3"), rob_var(x = x, y = y, type = "S3"))
  testthat::expect_equal(rob_var(x = y, y = x, type = "S4"), rob_var(x = x, y = y, type = "S4"))

  ## Check output ----

  # For 'type' = 'S1', the output should be 2
  testthat::expect_equal(rob_var(x = x, y = y, type = "S1"), 2)

  # For 'type' = 'S2', the output should be 2
  testthat::expect_equal(rob_var(x = x, y = y, type = "S2"), 2)

  # For 'type' = 'S3', the output should be 3
  testthat::expect_equal(rob_var(x = x, y = y, type = "S3"), 3)

  # For 'type' = 'S4', the output should be 3
  testthat::expect_equal(rob_var(x = x, y = y, type = "S4"), 3)

  # The output should be a numeric scalar
  checkmate::expect_numeric(rob_var(x = x, y = y, type = "S1"), len = 1)
  checkmate::expect_numeric(rob_var(x = x, y = y, type = "S2"), len = 1)
  checkmate::expect_numeric(rob_var(x = x, y = y, type = "S3"), len = 1)
  checkmate::expect_numeric(rob_var(x = x, y = y, type = "S4"), len = 1)
})


