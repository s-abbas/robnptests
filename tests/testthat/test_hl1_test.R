context("HL1 test")

testthat::test_that("hl1_test works correctly", {
  ## ___________________________________________________________________________
  ## Computation of test statistic
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(5)
  y <- rnorm(5)

  res.s1 <- as.numeric(hl1_test(x, y, method = "exact", scale = "S1")$statistic)
  res.s2 <- as.numeric(hl1_test(x, y, method = "exact", scale = "S2")$statistic)

  testthat::expect_equal(res.s1, (hodges_lehmann(x) - hodges_lehmann(y))/rob_var(x, y, type = "S1"))
  testthat::expect_equal(res.s2, (hodges_lehmann(x) - hodges_lehmann(y))/rob_var(x, y, type = "S2"))

  ## ___________________________________________________________________________
  ## If not method is given, the asymptotic test should be performed automatically
  ## for large sample sizes
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(30)
  y <- rnorm(30)

  estimates <- c(hodges_lehmann(x), hodges_lehmann(y - delta))

  m <- length(x)
  n <- length(y)

  # pairwise differences the density estimate is calculated from:
  xcomb <- utils::combn(x, 2)
  ycomb <- utils::combn(y - delta, 2)
  pwdiffs <- c(xcomb[2, ] - xcomb[1, ], ycomb[2, ] - ycomb[1, ])
  dens <- stats::density(pwdiffs)
  dens <- stats::approxfun(dens)

  int <- dens(0)

  res.asymptotic <- hl1_test(x, y)
  res.statistic <- as.numeric(res.asymptotic$statistic)
  res.method <- res.asymptotic$method

  testthat::expect_equal(res.statistic, sqrt(12*m*n/(m+n)) * int * (estimates[1] - estimates[2]))
  testthat::expect_equal(res.method, "Asymptotic test based on the Hodges-Lehmann estimator")

  ## ___________________________________________________________________________
  ## If not method is given, the asymptotic test should be performed automatically
  ## for small to moderate sample sizes
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(10)
  y <- rnorm(10)

  res.randomized <- hl1_test(x, y)
  res.method <- res.asymptotic$method

  testthat::expect_equal(res.method, "Randomization test based on the Hodges-Lehmann estimator")
}
)
