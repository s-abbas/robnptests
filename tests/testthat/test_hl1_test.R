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

  delta <- 0
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
  testthat::expect_equal(res.method, "Asymptotic test based on the one-sample Hodges-Lehmann estimator")

  ## ___________________________________________________________________________
  ## If no method is given, the randomized test should be performed automatically
  ## for small to moderate sample sizes
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(10)
  y <- rnorm(10)

  res.randomized <- hl1_test(x, y)
  res.method <- res.randomized$method

  testthat::expect_equal(res.method, "Randomization test based on the one-sample Hodges-Lehmann estimator")

  ## ___________________________________________________________________________
  ## Throw error when one or both samples consist of less than five observations
  ## ___________________________________________________________________________

  ## One sample is too small
  set.seed(108)

  x <- rnorm(4)
  y <- rnorm(20)

  testthat::expect_error(hl1_test(x, y))

  ## Many missing values
  x <- c(x, rep(NA, 5))

  testthat::expect_error(hl1_test(x, y, na.rm = TRUE))


  ## ___________________________________________________________________________
  ## Return errors if the wrong methods/alternative are handed over
  ## ___________________________________________________________________________

  set.seed(108)

  x <- rnorm(30)
  y <- rnorm(30)

  testthat::expect_error(hl1_test(x, y, alternative="none"))
  testthat::expect_error(hl1_test(x, y, method="toss_a_coin"))
  testthat::expect_error(hl1_test(x, y, delta="xyz"))
  testthat::expect_error(hl1_test(x, y, scale="SD"))

  ## ___________________________________________________________________________
  ## Compare p-values for different alternative hypothesis
  ## ___________________________________________________________________________
  ## The p-values should be related by the following equations:
  ## (i) p.two.sided = 2 * min(p.less, p.greater)
  ## (ii) p.less = 1 - p.greater,
  ## where p.two.sided is the p-value for the two.sided alternative and
  ## p.greater and p.less are the p-values for the one-sided alternatives.

  ##
  ## Asymptotic test
  ##

  set.seed(108)

  x <- rnorm(30)
  y <- rnorm(30)

  p.two.sided <- hl1_test(x, y, method = "asymptotic", alternative = "two.sided")$p.value
  p.greater <- hl1_test(x, y, method = "asymptotic", alternative = "greater")$p.value
  p.less <- hl1_test(x, y, method = "asymptotic", alternative = "less")$p.value

  testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater))
  testthat::expect_equal(p.less, 1 - p.greater)

  ##
  ## Permutation test
  ##

  set.seed(108)

  x <- rnorm(5)
  y <- rnorm(5)

  p.two.sided <- hl1_test(x, y, method = "exact", alternative = "two.sided")$p.value
  p.greater <- hl1_test(x, y, method = "exact", alternative = "greater")$p.value
  p.less <- hl1_test(x, y, method = "exact", alternative = "less")$p.value

  testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater))

  ##
  ## Randomized test
  ##

  set.seed(108)

  x <- rnorm(10)
  y <- rnorm(10)

  set.seed(123)
  p.two.sided <- hl1_test(x, y, method = "randomization", alternative = "two.sided", n.rep = 1000)$p.value
  set.seed(123)
  p.greater <- hl1_test(x, y, method = "randomization", alternative = "greater", n.rep = 1000)$p.value
  set.seed(123)
  p.less <- hl1_test(x, y, method = "randomization", alternative = "less", n.rep = 1000)$p.value

  ## We use the correction by Phipson and Smyth (2011) as implemented in the
  ## R package statmod. Thus, we need to determine the number b of values in
  ## randomization distributions which are at least as extreme as the observed
  ## value.
  set.seed(99)
  dist.sampled <- perm_distribution(x, y, type = "HL11", randomization = TRUE, n.rep = 1000)
  b <- sum(dist.sampled >= rob_perm_statistic(x, y, type = "HL11")$statistic)/1001

  testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater))

}
)
