context("HL2 test")

testthat::test_that("hl2_test works correctly", {
  ## ___________________________________________________________________________
  ## Computation of test statistic
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(5)
  y <- rnorm(5)

  res.s1 <- as.numeric(hl2_test(x, y, method = "permutation", scale = "S1")$statistic)
  res.s2 <- as.numeric(hl2_test(x, y, method = "permutation", scale = "S2")$statistic)

  testthat::expect_equal(res.s1, hodges_lehmann_2sample(x, y) / rob_var(x, y, type = "S1"))
  testthat::expect_equal(res.s2, hodges_lehmann_2sample(x, y) / rob_var(x, y, type = "S2"))


  ## ___________________________________________________________________________
  ## If not method is given, the asymptotic test should be performed automatically
  ## for large sample sizes
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(30)
  y <- rnorm(30)

  delta <- 0
  m <- length(x)
  n <- length(y)
  lambda <- m/(m + n)

  ## Estimation of density at zero for pairwise differences
  xcomb <- utils::combn(x, 2)
  ycomb <- utils::combn(y - delta, 2)

  pwdiffs <- c(xcomb[2, ] - xcomb[1, ], ycomb[2, ] - ycomb[1, ])

  dens <- stats::density(pwdiffs)
  int <- stats::approxfun(dens)(0)

  est <- hodges_lehmann_2sample(x, y - delta)
  statistic <- sqrt(12 * lambda * (1 - lambda)) * int * sqrt(m + n) * est

  # estimates <- hodges_lehmann_2sample(x, y)

  res.asymptotic <- hl2_test(x, y)
  res.statistic <- as.numeric(res.asymptotic$statistic)
  res.method <- res.asymptotic$method

  testthat::expect_equal(res.statistic, sqrt(12 * lambda * (1 - lambda)) * int * sqrt(m + n) * est)
  testthat::expect_equal(res.method, "Asymptotic test based on the Two-Sample Hodges-Lehmann estimator")

  ## ___________________________________________________________________________
  ## If no method is given, the randomized test should be performed automatically
  ## for small to moderate sample sizes
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(10)
  y <- rnorm(10)

  res.randomized <- hl2_test(x, y)
  res.method <- res.randomized$method

  testthat::expect_equal(res.method, "Randomization test based on the Two-Sample Hodges-Lehmann estimator")

  ## ___________________________________________________________________________
  ## Throw error when one or both samples consist of less than five observations
  ## ___________________________________________________________________________

  ## One sample is too small
  set.seed(108)

  x <- rnorm(4)
  y <- rnorm(20)

  testthat::expect_error(hl2_test(x, y))

  ## Many missing values
  x <- c(x, rep(NA, 5))

  testthat::expect_error(hl2_test(x, y, na.rm = TRUE))


  ## ___________________________________________________________________________
  ## Return errors if the wrong methods/alternative are handed over
  ## ___________________________________________________________________________

  set.seed(108)

  x <- rnorm(30)
  y <- rnorm(30)

  testthat::expect_error(hl2_test(x, y, alternative="none"))
  testthat::expect_error(hl2_test(x, y, method="toss_a_coin"))
  testthat::expect_error(hl2_test(x, y, delta="xyz"))
  testthat::expect_error(hl2_test(x, y, scale="SD"))

  ## ___________________________________________________________________________
  ## Compare p-values for different alternative hypothesis
  ## ___________________________________________________________________________
  ## The p-values should be related by the following equations:
  ## (i) p.two.sided = 2 * min(p.less, p.greater)
  ## (ii) p.less = 1 - p.greater,
  ## where p.two.sided is the p-value for the two.sided alternative and
  ## p.greater and p.less are the p-values for the one-sided alternatives.
  ##
  ## For the permutation and the randomization test, we need to increase the
  ## tolerance in the comparison. This is because the null distributions are
  ## discrete. Hence,

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

  p.two.sided <- hl2_test(x, y, method = "permutation", alternative = "two.sided")$p.value
  p.greater <- hl2_test(x, y, method = "permutation", alternative = "greater")$p.value
  p.less <- hl2_test(x, y, method = "permutation", alternative = "less")$p.value

  perm.dist <- perm_distribution(x, y, type = "HL21", randomization = FALSE)
  hl21.statistic <- rob_perm_statistic(x, y, type = "HL21")$statistic

  testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater))

  ## In the comparison of the one-sided p-values, we need to add the number of
  ## values in the permutation distribution that are equal to the value of the
  ## test statistic. Because of the discrete null distribution, the value of the
  ## test statistic is included in the computation of the left-sided and the
  ## computation of the right-sided p-value. Hence, it is counted twice so
  ## that p.less + p.greater > 1.
  testthat::expect_equal(p.less, 1 - p.greater + sum(hl21.statistic == perm.dist)/252)

  ##
  ## Randomization test
  ##

  set.seed(108)

  x <- rnorm(10)
  y <- rnorm(10)

  set.seed(123)
  p.two.sided <- hl2_test(x, y, method = "randomization", alternative = "two.sided", n.rep = 1000)$p.value
  set.seed(123)
  p.greater <- hl2_test(x, y, method = "randomization", alternative = "greater", n.rep = 1000)$p.value
  set.seed(123)
  p.less <- hl2_test(x, y, method = "randomization", alternative = "less", n.rep = 1000)$p.value

  ## We increase the tolerance for the comparisons. One reason is the discrete
  ## null distribution. Moreover, as we use the correction by Phipson and Smyth (2011),
  ## it would be necessary to compute and add the integral in equation (2)
  ## of their paper, which would make this test case more complicated.
  testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater), tolerance = 10^(-2))
  testthat::expect_equal(1, p.less + p.greater, tolerance = 10^(-3))

  ##
  ## Zeros in var.test: A warning should be thrown if at least 1 of the observations is zero.
  ##

  set.seed(108)

  x <- rnorm(10)
  y <- c(rnorm(9), 0)

  testthat::expect_warning(hl2_test(x, y, method = "asymptotic", var.test = TRUE))
}

)