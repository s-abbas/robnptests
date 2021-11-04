testthat::test_that("trimmed_test works correctly", {

  testthat::skip_on_cran()

  # Exemplary input vectors ----
  set.seed(108)

  x <- rnorm(30)
  y <- rnorm(30)

  # Create and compare snapshots of test output ----

  # Permutation test
  testthat::expect_snapshot_output(trimmed_test(x = x[1:5], y = y[1:5],
                                                method = "permutation"))

  testthat::expect_snapshot_output(trimmed_test(x = x[1:5], y = y[1:5],
                                                method = "permutation", var.test = TRUE))

  # Randomization test
  testthat::expect_snapshot_output(trimmed_test(x = x[1:10], y = y[1:10],
                                                method = "randomization", var.test = TRUE))

  testthat::expect_snapshot_output(trimmed_test(x = x[1:10], y = y[1:10],
                                                method = "randomization",
                                                n.rep = 10000, var.test = TRUE))

  # Asymptotic test
  testthat::expect_snapshot_output(trimmed_test(x = x, y = y, method = "asymptotic"))
  testthat::expect_snapshot_output(trimmed_test(x = x, y = y, method = "asymptotic",
                                                var.test = TRUE))

  # Compare value of the test statistic to manually computed value ----

  res <- as.numeric(trimmed_test(x = x, y = y, method = "asymptotic")$statistic)

  var.x <- win_var(x, gamma = 0.2)
  var.y <- win_var(y, gamma = 0.2)
  h.x <- var.x$h
  h.y <- var.y$h
  var.x <- var.x$var
  var.y <- var.y$var
  m <- length(x)
  n <- length(y)
  pool.var <- ((m - 1) * var.x + (n - 1) * var.y)/(h.x + h.y - 2)

  testthat::expect_equal(res, (trim_mean(x, gamma = 0.2) - trim_mean(y, gamma = 0.2))/sqrt(pool.var * (1/h.x + 1/h.y)))

  # Automatic selection of the method to compute the p-value ----

  # Asymptotic test for large samples
  testthat::expect_equal(trimmed_test(x = x, y = y)$method,
                         "Yuen's t-test")

  # Randomization test for small samples
  testthat::expect_equal(trimmed_test(x = x[1:10], y = y[1:10], n.rep = 100)$method,
                         "Randomization test based on trimmed means (100 random permutations)")

  # Permutation test if sample size is small and 'n.rep' equals the number of
  # possible splits
  testthat::expect_equal(trimmed_test(x = x[1:5], y = y[1:5], n.rep = 252)$method,
                         "Exact permutation test based on trimmed means")

  # User-specified selection of the method to compute the p-value ----

  # Asymptotic test
  testthat::expect_equal(trimmed_test(x = x, y = y, method = "asymptotic")$method,
                         "Yuen's t-test")

  # Randomization test for small samples
  testthat::expect_equal(trimmed_test(x = x[1:5], y = y[1:5], method = "randomization",
                                      n.rep = 100)$method,
                         "Randomization test based on trimmed means (100 random permutations)")

  # Permutation test if sample size is small and 'n.rep' equals the number of
  # possible splits
  testthat::expect_equal(trimmed_test(x = x[1:5], y = y[1:5], method = "permutation")$method,
                         "Exact permutation test based on trimmed means")

  # One of the sample contains less than five non-missing observations ----

  testthat::expect_error(trimmed_test(x = x[1:4], y = y))
  testthat::expect_error(trimmed_test(x = x, y = c(y[1:4], rep(NA_real_, 10))))

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
  p.two.sided <- trimmed_test(x = x, y = y, method = "asymptotic",
                              alternative = "two.sided")$p.value
  p.greater <- trimmed_test(x = x, y = y, method = "asymptotic",
                            alternative = "greater")$p.value
  p.less <- trimmed_test(x = x, y = y, method = "asymptotic",
                         alternative = "less")$p.value

  testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater))
  testthat::expect_equal(p.less, 1 - p.greater)

  # Permutation test
  p.two.sided <- trimmed_test(x = x[1:5], y = y[1:5], method = "permutation",
                              alternative = "two.sided")$p.value
  p.greater <- trimmed_test(x = x[1:5], y = y[1:5], method = "permutation",
                            alternative = "greater")$p.value
  p.less <- trimmed_test(x = x[1:5], y = y[1:5], method = "permutation",
                         alternative = "less")$p.value

  testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater))

  # In the comparison of the one-sided p-values, we need to add the number of
  # values in the permutation distribution that are equal to the value of the
  # test statistic. Because of the discrete null distribution, the value of the
  # test statistic is included in the computation of the left-sided and the
  # computation of the right-sided p-value. Hence, it is counted more than once
  # so that p.less + p.greater > 1.
  splits <- gtools::combinations(10, 5, 1:10)
  complete <- c(x[1:5], y[1:5])

  perm.dist <- apply(splits, 1, function(s) {
    trimmed_t(x = complete[s], y = complete[-s], gamma = 0.2)$statistic
  })

  trimmed.t.statistic <- trimmed_t(x = x[1:5], y = y[1:5])$statistic
  testthat::expect_equal(p.less, 1 - p.greater + sum(trimmed.t.statistic == perm.dist)/252)

  # Randomization test
  set.seed(168)
  p.two.sided <- trimmed_test(x = x[1:10], y = y[1:10], method = "randomization",
                              alternative = "two.sided", n.rep = 10000)$p.value
  set.seed(168)
  p.greater <- trimmed_test(x = x[1:10], y = y[1:10], method = "randomization",
                            alternative = "greater", n.rep = 10000)$p.value
  set.seed(168)
  p.less <- trimmed_test(x = x[1:10], y = y[1:10], method = "randomization",
                         alternative = "less", n.rep = 10000)$p.value

  # We increase the tolerance for the comparisons. One reason is the discrete
  # null distribution. Moreover, as we use the correction by Phipson and Smyth
  # (2011), it would be necessary to compute and add the integral in equation (2)
  # of their paper, which would make this test case more complicated.
  testthat::expect_true(abs(p.two.sided - 2 * min(p.less, p.greater)) < 10^(-2))
  testthat::expect_true(abs(1 - (p.less + p.greater)) < 10^(-2))

  # Test for scale difference ----

  # One of the samples contains zeros
  testthat::expect_message(trimmed_test(x = x[1:10], y = c(y[1:9], 0),
                                        method = "asymptotic", var.test = TRUE))
})
