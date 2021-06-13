# Computation of permutation distribution ----
testthat::test_that("perm_distribution works correctly", {

  # testthat::skip_on_cran()

  # Exemplary input vectors ----
  x <- 1:5
  y <- 6:10
  m <- length(x)
  n <- length(y)

  types <- c("HL11", "HL12", "HL21", "HL22", "MED1", "MED2")

  for (type in types) {

    ## Check output ----

    ## The output should be a numeric vector of the specified length for all tests

    # Permutation distribution
    checkmate::expect_numeric(perm_distribution(x = x, y = y, type = type), len = choose(m + n, m))

    # Check that it does not accidentally calculate a randomization distribution
    # if we e.g. hand over n.rep
    testthat::expect_equal(perm_distribution(x = x, y = y, type = type, randomization = FALSE),
                           perm_distribution(x = x, y = y, type = type, n.rep = 1000))
    # Randomization distribution
    checkmate::expect_numeric(perm_distribution(x = x, y = y, type = type, randomization = TRUE, n.rep = 100), len = 100)

    # Are randomization distributions reproducible?
    testthat::expect_equal(
      { set.seed(710); perm_distribution(x = x, y = y, type = type, randomization = TRUE, n.rep = 100) },
      { set.seed(710); perm_distribution(x = x, y = y, type = type, randomization = TRUE, n.rep = 100) }
    )
    testthat::expect_error(perm_distribution(x = x, y = y, type = type, randomization = TRUE, n.rep = 1000))
  }
})

# Computation of permutation distribution for M-test statistics ----
testthat::test_that("m_est_perm_distribution works correctly", {

  # testthat::skip_on_cran()

  # Exemplary input vectors ----
  x <- 1:5
  y <- 6:10
  m <- length(x)
  n <- length(y)

  psis <- c("huber", "hampel", "bisquare")
  ks <- list("huber" = 1.345, "hampel" = c(1, 2, 3), "bisquare" = 1.345)

  for (psi in psis) {

    # Check output ----

    # The output should be a numeric vector

    # Permutation distribution
    checkmate::expect_numeric(m_est_perm_distribution(x = x, y = y, psi = psi, k = ks[[psi]]), len = choose(m + n, m))
    # Randomization distribution
    checkmate::expect_numeric(m_est_perm_distribution(x = x, y = y, psi = psi, k = ks[[psi]], randomization = TRUE, n.rep = 100), len = 100)

    # Check that we actually get the permutation distribution if we expect it
    testthat::expect_equal(m_est_perm_distribution(x = x, y = y, psi = psi, k = ks[[psi]], randomization = FALSE),
                           m_est_perm_distribution(x = x, y = y, psi = psi, k = ks[[psi]], n.rep = 1000))

    # Check reproducibility of randomization distribution:
    testthat::expect_equal(
      { set.seed(710); m_est_perm_distribution(x = x, y = y, psi = psi, k = ks[[psi]], randomization = TRUE, n.rep = 100) },
      { set.seed(710); m_est_perm_distribution(x = x, y = y, psi = psi, k = ks[[psi]], randomization = TRUE, n.rep = 100) }
    )

    # Check that we get an error if n.rep exceeds choose(m + n, n)
    testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = psi, k = ks[[psi]], randomization = TRUE, n.rep = 1000))
  }

})

# Computation of permutation p-value ----
testthat::test_that("calc_perm_p_value works correctly", {

  # testthat::skip_on_cran()

  # Check output ----

  distribution <- 1:252
  statistic <- 50

  # The output should be a numeric scalar

  # Permutation distribution

  alternatives <- c("two.sided", "greater", "less")

  # Check that p-values are between 0 and 1
  for(alternative in alternatives) {

    checkmate::expect_number(calc_perm_p_value(
      statistic = statistic, distribution = distribution, m = 5, n = 5,
      randomization = FALSE, n.rep = 10000, alternative = alternative), lower = 0, upper = 1)

  }

  # Is the p-value we expect returned?
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
  # We trust the permp-function from statmod for calculation of the
  # randomization p-values
})
