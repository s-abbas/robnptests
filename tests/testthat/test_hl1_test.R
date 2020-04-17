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


  ##### THERE'S A DELTA MISSING HERE?

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
  ## If not method is given, the randomized test should be performed automatically
  ## for small to moderate sample sizes
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(10)
  y <- rnorm(10)

  res.randomized <- hl1_test(x, y)
  res.method <- res.randomized$method

  testthat::expect_equal(res.method, "Randomization test based on the Hodges-Lehmann estimator")


  ##
  ## Return errors if the wrong methods/alternative are handed over
  ##

  set.seed(108)

  x <- rnorm(30)
  y <- rnorm(30)

  testthat::expect_error(hl1_test(x, y, alternative="none"))
  testthat::expect_error(hl1_test(x, y, method="toss_a_coin"))
  testthat::expect_error(hl1_test(x, y, delta="xyz"))
  testthat::expect_error(hl1_test(x, y, scale="SD"))

  ##
  ## Results should differ depending on the alternative given to the test
  ##

  set.seed(108)

  x <- rnorm(30)
  y <- rnorm(30) + 0.5

  p.two.sided <- hl1_test(x, y, method="asymptotic", alternative="two.sided")$p.value
  p.greater <- hl1_test(x, y, method="asymptotic", alternative="greater")$p.value
  p.less <- hl1_test(x, y, method="asymptotic", alternative="less")$p.value

  testthat::expect_false(p.two.sided == p.greater)
  testthat::expect_false(p.greater == p.less)
  testthat::expect_false(p.two.sided == p.less)
}
)
