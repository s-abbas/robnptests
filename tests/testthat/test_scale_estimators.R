context("Scale Estimators")

testthat::test_that("Robust scale estimators work correctly", {

  ##
  ## Check that the function returns errors when we hand over the wrong method
  ## or something non-numeric
  ##

  set.seed(108)

  x <- rnorm(10)
  y <- rnorm(10)

  testthat::expect_error(rob_var(x, y, na.rm=FALSE, type="not implemented"))
  testthat::expect_error(rob_var(x, y, na.rm=FALSE, type=1))
  testthat::expect_error(rob_var(letters[1:10], y, na.rm=FALSE))

  ##
  ## NA handling
  ##

  types <- c("S1", "S2", "S3", "S4")

  x2 <- c(x, NA)

  # invisible(sapply(types, function(type) {
  #   testthat::expect_true(is.na(rob_var(x2, y, na.rm=FALSE, type=type)))
  #   testthat::expect_identical(rob_var(x, y, na.rm=FALSE, type=type),
  #                              rob_var(x2, y, na.rm=TRUE, type=type))
  # })) # Kann man das so machen, um alle vier Schätzer gleichzeitig abzuarbeiten?
      # zumindest gibt es einen error aus wenn sie nicht identical sind.
      # invisible sorgt dafür, dass kein Output in die Konsole geprintet wird

  ##
  ## Check for location and permutation invariance
  ##

  x3 <- x + 5 # shift in location
  x4 <- sort(x) # a simple 'permutation' of x

  # invisible(sapply(types, function(type) {
  #   testthat::expect_equal(rob_var(x, y, type=type), rob_var(y, x, type=type)) # switching x and y
  #   testthat::expect_equal(rob_var(x, y, type=type), rob_var(x3, y, type=type)) # location invariance
  #   testthat::expect_equal(rob_var(x, y, type=type), rob_var(x4, y, type=type)) # permutation invariance
  # }))

  ##
  ## Check if the function returns 0 when all observations are constant
  ## and prints a warning if the observations are not.
  ##

  # invisible(sapply(types, function(type) {-
  #   testthat::expect_equal(suppressWarnings(rob_var(rep(1, 10), rep(1, 10), type=type)), 0)
  #   testthat::expect_warning(rob_var(x, rep(1, 10), type=type))
  # }))

  testthat::expect_error(rob_var(rep(c(1, 2), 5), c(rep(1, 9), 100), type="S1"))

}
)


testthat::test_that("Winsorized variance works correctly", {

  ##
  ## Check that the function returns errors when we hand over the wrong parameter
  ## or something non-numeric
  ##

  set.seed(108)

  x <- rnorm(10)
  y <- rnorm(10)

  testthat::expect_error(win_var(x, gamma=-0.2))
  testthat::expect_error(win_var(x, gamma="none"))

  ##
  ## NA handling
  ##

  x2 <- c(x, NA)

  testthat::expect_true(is.na(win_var(x2, gamma=0)))
  testthat::expect_equal(win_var(x, gamma=0), win_var(x2, gamma=0, na.rm=TRUE))
  testthat::expect_equal(win_var(x, gamma=0.25), win_var(x2, gamma=0.25, na.rm=TRUE))

  ##
  ## Check for location and permutation invariance
  ##

  x3 <- x + 5
  x4 <- sort(x)

  testthat::expect_equal(win_var(x, gamma=0), win_var(x3, gamma=0, na.rm=TRUE))
  testthat::expect_equal(win_var(x, gamma=0), win_var(x4, gamma=0, na.rm=TRUE))
  testthat::expect_equal(win_var(x, gamma=0.25), win_var(x3, gamma=0.25, na.rm=TRUE))
  testthat::expect_equal(win_var(x, gamma=0.25), win_var(x4, gamma=0.25, na.rm=TRUE))

  ##
  ## Check whether different values of gamma yield different estimates
  ##

  g1 <- win_var(x, gamma=0)$var
  g2 <- win_var(x, gamma=0.1)$var

  testthat::expect_true(g1 != g2)

})

# testthat::test_that("Asymmetric trimmed variance works correctly", {
#
#  ## TO-DO
#
# })