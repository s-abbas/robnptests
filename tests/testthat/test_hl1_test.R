testthat::test_that("hl1_test works correctly", {

  # testthat::skip_on_cran()
  # Exemplary input vectors ----
  set.seed(108)

  x <- rnorm(30)
  y <- rnorm(30)

  # Create and compare snapshots of test output ----

  # Permutation test
  testthat::expect_snapshot_output(hl1_test(x = x[1:5], y = y[1:5],
                                            method = "permutation", scale = "S1"))
  testthat::expect_snapshot_output(hl1_test(x = x[1:5], y = y[1:5],
                                            method = "permutation", scale = "S2"))

  testthat::expect_snapshot_output(hl1_test(x = x[1:5], y = y[1:5],
                                            method = "permutation", scale = "S1",
                                            var.test = TRUE))
  testthat::expect_snapshot_output(hl1_test(x = x[1:5], y = y[1:5],
                                            method = "permutation", scale = "S2",
                                            var.test = TRUE))

  # Randomization test
  testthat::expect_snapshot_output(hl1_test(x = x[1:10], y = y[1:10],
                                            method = "randomization",
                                            n.rep =  10000, scale = "S1"))
  testthat::expect_snapshot_output(hl1_test(x = x[1:10], y = y[1:10],
                                            method = "randomization",
                                            n.rep = 10000,  scale = "S2"))

  testthat::expect_snapshot_output(hl1_test(x = x[1:10], y = y[1:10],
                                            method = "randomization",
                                            n.rep = 10000,  scale = "S1",
                                            var.test = TRUE))
  testthat::expect_snapshot_output(hl1_test(x = x[1:10], y = y[1:10],
                                            method = "randomization",
                                            n.rep = 10000,  scale = "S2",
                                            var.test = TRUE))

  # Asymptotic test
  testthat::expect_snapshot_output(hl1_test(x = x, y = y, method = "asymptotic"))
  testthat::expect_snapshot_output(hl1_test(x = x, y = y, method = "asymptotic",
                                            var.test = TRUE))

  # Compare value of the test statistic to manually computed value ----

  res.s1 <- as.numeric(hl1_test(x = x, y = y, method = "randomization", n.rep = 100,
                                scale = "S1")$statistic)
  res.s2 <- as.numeric(hl1_test(x = x, y = y, method = "randomization", n.rep = 100,
                                scale = "S2")$statistic)

  testthat::expect_equal(res.s1,
                         (hodges_lehmann(x) - hodges_lehmann(y))/rob_var(x, y, type = "S1"))
  testthat::expect_equal(res.s2,
                         (hodges_lehmann(x) - hodges_lehmann(y))/rob_var(x, y, type = "S2"))

  # Automatic selection of the method to compute the p-value ----

  # Asymptotic test for large samples
  testthat::expect_equal(hl1_test(x = x, y = y)$method,
                         "Asymptotic test based on HL1-estimator")

  # Randomization test for small samples
  testthat::expect_equal(hl1_test(x = x[1:10], y = y[1:10], n.rep = 100)$method,
                         "Randomization test based on HL1-estimator (100 random permutations)")

  # Permutation test if sample size is small and 'n.rep' equals the number of
  # possible splits
  testthat::expect_equal(hl1_test(x = x[1:5], y = y[1:5], n.rep = 252)$method,
                         "Exact permutation test based on HL1-estimator")

  # User-specified selection of the method to compute the p-value ----

  # Asymptotic test
  testthat::expect_equal(hl1_test(x = x, y = y, method = "asymptotic")$method,
                         "Asymptotic test based on HL1-estimator")

  # Randomization test for small samples
  testthat::expect_equal(hl1_test(x = x[1:5], y = y[1:5], method = "randomization",
                                  n.rep = 100)$method,
                         "Randomization test based on HL1-estimator (100 random permutations)")

  # Permutation test for small samples
  testthat::expect_equal(hl1_test(x = x[1:5], y = y[1:5], method = "permutation")$method,
                         "Exact permutation test based on HL1-estimator")

  # One of the sample contains less than five non-missing observations ----

  testthat::expect_error(hl1_test(x = x[1:4], y = y))
  testthat::expect_error(hl1_test(x = x, y = c(y[1:4], rep(NA_real_, 10))))

  # Computation of the p-values ----
  # The p-values should be related by the following equations:
  # (i) p.two.sided = 2 * min(p.less, p.greater)
  # (ii) p.less = 1 - p.greater,
  # where p.two.sided is the p-value for the two.sided alternative and
  # p.greater and p.less are the p-values for the one-sided alternatives.
  #
  # For the permutation and the randomization test, we need to increase the
  # tolerance in the comparison. This is because the null distributions are
  # discrete.

  # Asymptotic test
  p.two.sided <- hl1_test(x = x, y = y, method = "asymptotic",
                          alternative = "two.sided")$p.value
  p.greater <- hl1_test(x = x, y = y, method = "asymptotic",
                        alternative = "greater")$p.value
  p.less <- hl1_test(x = x, y = y, method = "asymptotic",
                     alternative = "less")$p.value

  testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater))
  testthat::expect_equal(p.less, 1 - p.greater)

  # Permutation test
  p.two.sided <- hl1_test(x = x[1:5], y = y[1:5], method = "permutation",
                          alternative = "two.sided")$p.value
  p.greater <- hl1_test(x = x[1:5], y = y[1:5], method = "permutation",
                        alternative = "greater")$p.value
  p.less <- hl1_test(x = x[1:5], y = y[1:5], method = "permutation",
                     alternative = "less")$p.value

  testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater))

  # In the comparison of the one-sided p-values, we need to add the number of
  # values in the permutation distribution that are equal to the value of the
  # test statistic. Because of the discrete null distribution, the value of the
  # test statistic is included in the computation of the left-sided and the
  # computation of the right-sided p-value. Hence, it is counted than once
  # so that p.less + p.greater > 1
  perm.dist <- perm_distribution(x = x[1:5], y = y[1:5], type = "HL11", randomization = FALSE)
  hl11.statistic <- rob_perm_statistic(x = x[1:5], y = y[1:5], type = "HL11")$statistic
  testthat::expect_equal(p.less, 1 - p.greater + sum(hl11.statistic == perm.dist)/252)

  # Randomization test
  set.seed(168)
  p.two.sided <- hl1_test(x = x[1:10], y = y[1:10], method = "randomization",
                          alternative = "two.sided", n.rep = 10000)$p.value
  set.seed(168)
  p.greater <- hl1_test(x = x[1:10], y = y[1:10], method = "randomization",
                        alternative = "greater", n.rep = 10000)$p.value
  set.seed(168)
  p.less <- hl1_test(x = x[1:10], y = y[1:10], method = "randomization",
                     alternative = "less", n.rep = 10000)$p.value

  # We increase the tolerance for the comparisons. One reason is the discrete
  # null distribution. Moreover, as we use the correction by Phipson and Smyth
  # (2011), it would be necessary to compute and add the integral in equation (2)
  # of their paper, which would make this test case more complicated.
  testthat::expect_true(abs(p.two.sided - 2 * min(p.less, p.greater)) < 10^(-2))
  testthat::expect_true(abs(1 - (p.less + p.greater)) < 10^(-2))

  # Test for scale difference ----

  # One of the samples contains zeros
  testthat::expect_warning(hl1_test(x = x[1:10], y = c(y[1:9], 0),
                                    method = "asymptotic", var.test = TRUE))

  # Wobbling ----

  # Exemplary input vectors
  x <- c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
  y <- c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0)

  # Default is 'wobble' = FALSE so that an error is thrown
  testthat::expect_error(
    suppressWarnings(
      hl1_test(x, y, method = "randomization", n.rep = 1000)
    )
  )

  # Setting 'wobble' = TRUE only causes a warning
  testthat::expect_warning(hl1_test(x = x, y = y, method = "randomization",
                                    n.rep = 1000, wobble = TRUE,
                                    wobble.seed = 1234))

  # Comparison of wobbled values by 'hl1_test' and output of function 'wobble'
  set.seed(1234)
  wob <- wobble(x = x, y = y, check = FALSE)

  testthat::expect_equal(suppressWarnings(hl1_test(x = x, y = y,
                                                   method = "randomization",
                                                   n.rep = 1000,
                                                   wobble = TRUE,
                                                   wobble.seed = 1234)$statistic),
                         suppressWarnings(hl1_test(x = wob$x, y = wob$y,
                                                   method = "randomization",
                                                   n.rep = 1000,
                                                   wobble = FALSE)$statistic))
})
