
## Computation of permutation distribution ----
testthat::test_that("perm_distribution works correctly", {

  # testthat::skip_on_cran()

  ## Generate exemplary input vectors
  x <- 1:5
  y <- 6:10
  m <- length(x)
  n <- length(y)

  ## Check output ----

  ## The output should be a numeric vector of the specified length for all tests

  # Permutation distribution
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL11"), len = choose(m + n, m))
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL12"), len = choose(m + n, m))
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL21"), len = choose(m + n, m))
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL22"), len = choose(m + n, m))
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "MED1"), len = choose(m + n, m))
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "MED2"), len = choose(m + n, m))

  # Check that it does not accidentally calculate a randomization distribution if
  # we e.g. hand over n.rep
  testthat::expect_equal(perm_distribution(x = x, y = y, type = "HL11", randomization = FALSE),
                         perm_distribution(x = x, y = y, type = "HL11", n.rep = 1000))
  testthat::expect_equal(perm_distribution(x = x, y = y, type = "HL12", randomization = FALSE),
                         perm_distribution(x = x, y = y, type = "HL12", n.rep = 1000))
  testthat::expect_equal(perm_distribution(x = x, y = y, type = "HL21", randomization = FALSE),
                         perm_distribution(x = x, y = y, type = "HL21", n.rep = 1000))
  testthat::expect_equal(perm_distribution(x = x, y = y, type = "HL22", randomization = FALSE),
                         perm_distribution(x = x, y = y, type = "HL22", n.rep = 1000))
  testthat::expect_equal(perm_distribution(x = x, y = y, type = "MED1", randomization = FALSE),
                         perm_distribution(x = x, y = y, type = "MED1", n.rep = 1000))
  testthat::expect_equal(perm_distribution(x = x, y = y, type = "MED2", randomization = FALSE),
                         perm_distribution(x = x, y = y, type = "MED2", n.rep = 1000))

  # Randomization distribution
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL11", randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL12", randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL21", randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL22", randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "MED1", randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "MED2", randomization = TRUE, n.rep = 100), len = 100)

  # Are randomization distributions reproducible? I.e. do we get the same randomization distribution
  # if we set the same seed?
  testthat::expect_equal(
    { set.seed(710); perm_distribution(x = x, y = y, type = "HL11", randomization = TRUE, n.rep = 100) },
    { set.seed(710); perm_distribution(x = x, y = y, type = "HL11", randomization = TRUE, n.rep = 100) }
  )
  testthat::expect_equal(
    { set.seed(710); perm_distribution(x = x, y = y, type = "HL12", randomization = TRUE, n.rep = 100) },
    { set.seed(710); perm_distribution(x = x, y = y, type = "HL12", randomization = TRUE, n.rep = 100) }
  )
  testthat::expect_equal(
    { set.seed(710); perm_distribution(x = x, y = y, type = "HL21", randomization = TRUE, n.rep = 100) },
    { set.seed(710); perm_distribution(x = x, y = y, type = "HL21", randomization = TRUE, n.rep = 100) }
  )
  testthat::expect_equal(
    { set.seed(710); perm_distribution(x = x, y = y, type = "HL22", randomization = TRUE, n.rep = 100) },
    { set.seed(710); perm_distribution(x = x, y = y, type = "HL22", randomization = TRUE, n.rep = 100) }
  )
  testthat::expect_equal(
    { set.seed(710); perm_distribution(x = x, y = y, type = "MED1", randomization = TRUE, n.rep = 100) },
    { set.seed(710); perm_distribution(x = x, y = y, type = "MED1", randomization = TRUE, n.rep = 100) }
  )
  testthat::expect_equal(
    { set.seed(710); perm_distribution(x = x, y = y, type = "MED2", randomization = TRUE, n.rep = 100) },
    { set.seed(710); perm_distribution(x = x, y = y, type = "MED2", randomization = TRUE, n.rep = 100) }
  )
})

## Computation of permutation distribution for M-test statistics ----
testthat::test_that("m_est_perm_distribution works correctly", {

  # testthat::skip_on_cran()

  ## Generate exemplary input vectors
  x <- 1:5
  y <- 6:10
  m <- length(x)
  n <- length(y)

  ## Check output ----

  ## The output should be a numeric vector

  # Permutation distribution
  checkmate::expect_numeric(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345), len = choose(m + n, m))
  checkmate::expect_numeric(m_est_perm_distribution(x = x, y = y, psi = "hampel", k = c(1, 2, 3)), len = choose(m + n, m))
  checkmate::expect_numeric(m_est_perm_distribution(x = x, y = y, psi = "bisquare", k = 1.345), len = choose(m + n, m))

  # Randomization distribution
  checkmate::expect_numeric(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(m_est_perm_distribution(x = x, y = y, psi = "hampel", k = c(1, 2, 3), randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(m_est_perm_distribution(x = x, y = y, psi = "bisquare", k = 1.345, randomization = TRUE, n.rep = 100), len = 100)
})

## Computation of permutation p-value ----
testthat::test_that("calc_perm_p_value works correctly", {

  # testthat::skip_on_cran()

  ## Check output ----

  distribution <- 1:252
  statistic <- 50

  ## The output should be a numeric scalar

  # Permutation distribution

  # Check that p-values are between 0 and 1
  checkmate::expect_number(calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5, n = 5,
    randomization = FALSE, n.rep = 10000, alternative = "two.sided"), lower = 0, upper = 1)

  checkmate::expect_number(calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5, n = 5,
    randomization = FALSE, n.rep = 10000, alternative = "greater"), lower = 0, upper = 1)

  checkmate::expect_number(calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5, n = 5,
    randomization = FALSE, n.rep = 10000, alternative = "less"), lower = 0, upper = 1)

  # Is the p-value we expect returned
  testthat::expect_equal(calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5,
    n = 5, randomization = FALSE, n.rep = 100, alternative = "less"), 50 / 252)

  testthat::expect_equal(calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5,
    n = 5, randomization = FALSE, n.rep = 100, alternative = "greater"), 203 / 252)

  testthat::expect_equal(calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5,
    n = 5, randomization = FALSE, n.rep = 100, alternative = "two.sided"), 203 / 252)


  # Do p.greater and p.less sum up to 1?
  testthat::expect_equal(calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5,
    n = 5, randomization = FALSE, n.rep = 100, alternative = "less") + calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5,
    n = 5, randomization = FALSE, n.rep = 100, alternative = "greater") - 1/252, 1)

  # Does the relationship p.two.sided = 2 * min(p.greater, p.less) hold?
  # In this case it doesn't as the permutation distribution is not symmetric

  # Test with a symmetric distribution. Assume n = 2 and m = 8 (in this case the
  # number of splits is (10 over 2 = 45)

  distribution2 <- seq(-22, 22, 1)
  statistic2 <- 5

  p.less <- calc_perm_p_value(
    statistic = statistic2, distribution = distribution2, m = 2, n = 8,
    randomization = FALSE, n.rep = 100, alternative = "less")

  p.greater <- calc_perm_p_value(
    statistic = statistic2, distribution = distribution2, m = 2, n = 8,
    randomization = FALSE, n.rep = 100, alternative = "greater")

  p.two.sided <- calc_perm_p_value(
    statistic = statistic2, distribution = distribution2, m = 2, n = 8,
    randomization = FALSE, n.rep = 100, alternative = "two.sided")

  testthat::expect_equal(2 * min(p.less, p.greater), p.two.sided)


  # Randomization distribution
  # We trust the permp-function from statmod for calculation of the randomization p-value

  # Check that the output is between 0 and 1
  checkmate::expect_number(calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5, n = 5,
    randomization = TRUE, n.rep = 250, alternative = "two.sided"), lower = 0, upper = 1)
  checkmate::expect_number(calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5, n = 5,
    randomization = TRUE, n.rep = 250, alternative = "greater"), lower = 0, upper = 1)
  checkmate::expect_number(calc_perm_p_value(
    statistic = statistic, distribution = distribution, m = 5, n = 5,
    randomization = TRUE, n.rep = 250, alternative = "less"), lower = 0, upper = 1)
})


