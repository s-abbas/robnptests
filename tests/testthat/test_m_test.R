testthat::test_that("m_test works correctly", {
  #testthat::skip_on_cran()

  psi.funs <- c("huber", "hampel", "bisquare")

  # Exemplary input vectors ----
  set.seed(108)

  x <- rnorm(30)
  y <- rnorm(30)

  for (i in seq_along(psi.funs)) {
    psi <- psi.funs[i]

    # Create and compare snapshots of test output ----

    # Permutation test
    testthat::expect_snapshot_output(m_test(x = x[1:5], y = y[1:5],
                                            method = "permutation",
                                            psi = psi))

    testthat::expect_snapshot_output(m_test(x = x[1:5], y = y[1:5],
                                            method = "permutation",
                                            psi = psi,
                                            var.test = TRUE))

    # Randomization test
    testthat::expect_snapshot_output(m_test(x = x[1:10], y = y[1:10],
                                            psi = psi,
                                            method = "randomization",
                                            n.rep = 10000))

    testthat::expect_snapshot_output(m_test(x = x[1:10], y = y[1:10],
                                            psi = psi,
                                            method = "randomization",
                                            n.rep = 10000, var.test = TRUE))

    # Asymptotic test
    testthat::expect_snapshot_output(m_test(x = x, y = y,
                                            method = "asymptotic",
                                            psi = psi))
    testthat::expect_snapshot_output(m_test(x = x, y = y,
                                            method = "asymptotic",
                                            psi = psi,
                                            var.test = TRUE))

    # Automatic selection of the method to compute the p-value ----

    # Asymptotic test for large samples
    testthat::expect_equal(m_test(x = x, y = y, psi = psi)$method,
                           paste("Asymptotic test based on", paste0(toupper(substring(psi, 1, 1)), substring(psi, 2, nchar(psi))), "M-estimator"))

    # Randomization test for small samples
    testthat::expect_equal(m_test(x = x[1:10], y = y[1:10], psi = psi, n.rep = 100)$method,
                           paste("Randomization test based on", paste0(toupper(substring(psi, 1, 1)), substring(psi, 2, nchar(psi))), "M-estimator", paste0("(", 100), "random permutations)"))

    # Permutation test if sample size is small and 'n.rep' equals the number of
    # possible splits
    testthat::expect_equal(m_test(x = x[1:5], y = y[1:5], psi = psi, n.rep = 252)$method,
                           paste("Exact permutation test based on", paste0(toupper(substring(psi, 1, 1)), substring(psi, 2, nchar(psi))), "M-estimator"))

    # User-specified selection of the method to compute the p-value ----

    # Asymptotic test
    testthat::expect_equal(m_test(x = x, y = y, method = "asymptotic", psi = psi)$method,
                           paste("Asymptotic test based on", paste0(toupper(substring(psi, 1, 1)), substring(psi, 2, nchar(psi))), "M-estimator"))

    # Randomization test for small samples
    testthat::expect_equal(m_test(x = x[1:5], y = y[1:5], method = "randomization",
                                  psi = psi, n.rep = 100)$method,
                           paste("Randomization test based on", paste0(toupper(substring(psi, 1, 1)), substring(psi, 2, nchar(psi))), "M-estimator", paste0("(", 100), "random permutations)"))

    # Permutation test if sample size is small and 'n.rep' equals the number of
    # possible splits
    testthat::expect_equal(m_test(x = x[1:5], y = y[1:5], method = "permutation", psi = psi)$method,
                           paste("Exact permutation test based on", paste0(toupper(substring(psi, 1, 1)), substring(psi, 2, nchar(psi))), "M-estimator"))

    # One of the sample contains less than five non-missing observations ----

    testthat::expect_error(m_test(x = x[1:4], y = y, psi = psi))
    testthat::expect_error(
      suppressWarnings(
        m_test(x = x, y = c(y[1:4], rep(NA_real_, 10)), psi = psi)
      )
    )

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
    p.two.sided <- m_test(x = x, y = y, method = "asymptotic", psi = psi,
                          alternative = "two.sided")$p.value
    p.greater <- m_test(x = x, y = y, method = "asymptotic", psi = psi,
                        alternative = "greater")$p.value
    p.less <- m_test(x = x, y = y, method = "asymptotic", psi = psi,
                     alternative = "less")$p.value

    testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater))
    testthat::expect_equal(p.less, 1 - p.greater)

    # Permutation test
    p.two.sided <- m_test(x = x[1:5], y = y[1:5], method = "permutation", psi = psi,
                          alternative = "two.sided")$p.value
    p.greater <- m_test(x = x[1:5], y = y[1:5], method = "permutation", psi = psi,
                        alternative = "greater")$p.value
    p.less <- m_test(x = x[1:5], y = y[1:5], method = "permutation", psi = psi,
                     alternative = "less")$p.value

    testthat::expect_equal(p.two.sided, 2 * min(p.less, p.greater))

    # In the comparison of the one-sided p-values, we need to add the number of
    # values in the permutation distribution that are equal to the value of the
    # test statistic. Because of the discrete null distribution, the value of the
    # test statistic is included in the computation of the left-sided and the
    # computation of the right-sided p-value. Hence, it is counted more than once
    # so that p.less + p.greater > 1.
    perm.dist <- m_est_perm_distribution(x = x[1:5], y = y[1:5], randomization = FALSE,
                                         psi = psi, k = robustbase::.Mpsi.tuning.default(psi))

    m.statistic <- m_test_statistic(x = x[1:5], y = y[1:5], psi = psi)$statistic
    testthat::expect_equal(p.less, 1 - p.greater + sum(m.statistic == perm.dist)/252)

    # Randomization test
    set.seed(168)
    p.two.sided <- m_test(x = x[1:10], y = y[1:10], method = "randomization",
                          psi = psi, alternative = "two.sided", n.rep = 10000)$p.value
    set.seed(168)
    p.greater <- m_test(x = x[1:10], y = y[1:10], method = "randomization", psi = psi,
                        alternative = "greater", n.rep = 10000)$p.value
    set.seed(168)
    p.less <- m_test(x = x[1:10], y = y[1:10], method = "randomization", psi = psi,
                     alternative = "less", n.rep = 10000)$p.value

    # We increase the tolerance for the comparisons. One reason is the discrete
    # null distribution. Moreover, as we use the correction by Phipson and Smyth
    # (2011), it would be necessary to compute and add the integral in equation (2)
    # of their paper, which would make this test case more complicated.
    testthat::expect_true(abs(p.two.sided - 2 * min(p.less, p.greater)) < 10^(-2))
    testthat::expect_true(abs(1 - (p.less + p.greater)) < 10^(-2))

    # Test for scale difference ----

    # One of the samples contains zeros
    testthat::expect_message(m_test(x = x[1:10], y = c(y[1:9], 0),
                                    method = "asymptotic", psi = psi, var.test = TRUE))
  }
})
