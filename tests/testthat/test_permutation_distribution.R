
## Computation of permutation distribution ----
testthat::test_that("perm_distribution works correctly", {

  # testthat::skip_on_cran()

  ## Generate exemplary input vectors
  x <- 1:5
  y <- 6:10
  m <- length(x)
  n <- length(y)

  ## Check output ----

  ## The output should be a numeric vector

  # Permutation distribution
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL11"), len = choose(m + n, m))
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL12"), len = choose(m + n, m))
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL21"), len = choose(m + n, m))
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL22"), len = choose(m + n, m))
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "MED1"), len = choose(m + n, m))
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "MED2"), len = choose(m + n, m))

  # Randomization distribution
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL11", randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL12", randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL21", randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "HL22", randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "MED1", randomization = TRUE, n.rep = 100), len = 100)
  checkmate::expect_numeric(perm_distribution(x = x, y = y, type = "MED2", randomization = TRUE, n.rep = 100), len = 100)
})

## Computation of permutation distribution for M-test statistics ----
testthat::test_that("m_est_perm_distribution works correctly", {

  # testthat::skip_on_cran()

  ## Generate exemplary input vectors
  x <- 1:5
  y <- 6:10
  m <- length(x)
  n <- length(y)

  ## Checks for input arguments ----
  # 'x'
  testthat::expect_error(m_est_perm_distribution(y = y, psi = "huber", k = 1.345),
                         regexp = "'x' is missing.")
  testthat::expect_error(m_est_perm_distribution(x = c(x, NA), y = y, psi = "huber", k = 1.345),
                         regexp = "Assertion on 'x' failed: Contains missing values.")
  testthat::expect_error(m_est_perm_distribution(x = NULL, y = y, psi = "huber", k = 1.345),
                         regexp = "Assertion on 'x' failed: Must be of type 'numeric', not 'NULL'.")
  testthat::expect_error(m_est_perm_distribution(x = c(x, Inf), y = y, psi = "huber", k = 1.345),
                         regexp = "Assertion on 'x' failed: Must be finite.")
  testthat::expect_error(m_est_perm_distribution(x = x[1:2], y = y, psi = "huber", k = 1.345),
                         regexp = "Assertion on 'x' failed: Must have length >= 5, but has length 2.")

  # 'y'
  testthat::expect_error(m_est_perm_distribution(x = x, psi = "huber", k = 1.345),
                         regexp = "'y' is missing.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = c(y, NA), psi = "huber", k = 1.345),
                         regexp = "Assertion on 'y' failed: Contains missing values.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = NULL, psi = "huber", k = 1.345),
                         regexp = "Assertion on 'y' failed: Must be of type 'numeric', not 'NULL'.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = c(-Inf, y), psi = "huber", k = 1.345),
                         regexp = "Assertion on 'y' failed: Must be finite.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = y[3:5], psi = "huber", k = 1.345),
                         regexp = "Assertion on 'y' failed: Must have length >= 5, but has length 3.")

  # 'psi'
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, k = 1.345),
                         regexp = "'psi' is missing.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = 1, k = 1.345),
                         regexp = "Assertion on 'psi' failed: Must be element of set {'huber','hampel','bisquare'}, but types do not match (numeric != character).", fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = c("huber", "hampel"), k = 1.345),
                         regexp = "Assertion on 'psi' failed: Must be element of set {'huber','hampel','bisquare'}, but is not atomic scalar.", fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "tukey", k = 1.345), regexp = "Assertion on 'psi' failed: Must be element of set {'huber','hampel','bisquare'}, but is 'tukey'.", fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = NA, k = 1.345), regexp = "Assertion on 'psi' failed: Must be element of set {'huber','hampel','bisquare'}, but is 'NA'.", fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = NULL, k = 1.345), regexp = "Assertion on 'psi' failed: Must be a subset of {'huber','hampel','bisquare'}, not 'NULL'.", fixed = TRUE)

  # 'k'
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber"),
                         regexp = "'k' is missing.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = NA),
                         regexp = "Assertion on 'k' failed: Contains missing values.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = NULL),
                         regexp = "Assertion on 'k' failed: Must be of type 'numeric', not 'NULL'.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = Inf),
                         regexp = "Assertion on 'k' failed: Must be finite.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = c(1, 1.345)),
                         regexp = "Assertion on 'k' failed: Must have length 1, but has length 2.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "hampel", k = 1.345),
                         regexp = "Assertion on 'k' failed: Must have length 3, but has length 1.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = -1),
                         regexp = "Assertion on 'k' failed: Element 1 is not >= 0.")
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "hampel", k = c(1, 1.2, -1)),
                         regexp = "Assertion on 'k' failed: Element 3 is not >= 0.")

  # 'randomization'
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = NA),
                         regexp = "Assertion on 'randomization' failed: May not be NA.",
                         fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = 1),
                         regexp = "Assertion on 'randomization' failed: Must be of type 'logical flag', not 'double'.",
                         fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = c(TRUE, TRUE)),
                         regexp = "Assertion on 'randomization' failed: Must have length 1.",
                         fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = NULL),
                         regexp = "Assertion on 'randomization' failed: Must be of type 'logical flag', not 'NULL'.",
                         fixed = TRUE)

  # 'n.rep'
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = FALSE, n.rep = NA),
                         regexp = "Assertion on 'n.rep' failed: May not be NA.",
                         fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = FALSE, n.rep = "10000"),
                         regexp = "Assertion on 'n.rep' failed: Must be of type 'count', not 'character'.",
                         fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = FALSE, n.rep = c(1, 100)),
                         regexp = "Assertion on 'n.rep' failed: Must have length 1.",
                         fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = FALSE, n.rep = NULL),
                         regexp = "Assertion on 'n.rep' failed: Must be of type 'count', not 'NULL'.",
                         fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = FALSE, n.rep = -100),
                         regexp = "Assertion on 'n.rep' failed: Must be >= 1.",
                         fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = FALSE, n.rep = Inf),
                         regexp = "Assertion on 'n.rep' failed: Must be of type 'count', not 'double'.",
                         fixed = TRUE)
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345, randomization = TRUE, n.rep = 1000),
                         regexp = "'n.rep' may not be larger than 252, the number of all splits.",
                         fixed = TRUE)


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

  ## Generate exemplary input vectors
  statistic <- 50
  distribution <- 1:252

  ## Checks for input arguments ----
  # 'statistic'
  testthat::expect_error(calc_perm_p_value(distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "'statistic' is missing.")
  testthat::expect_error(calc_perm_p_value(statistic = "1", distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'statistic' failed: Must be of type 'number', not 'character'.")
  testthat::expect_error(calc_perm_p_value(statistic = NA, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'statistic' failed: May not be NA.")
  testthat::expect_error(calc_perm_p_value(statistic = NULL, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'statistic' failed: Must be of type 'number', not 'NULL'.")
  testthat::expect_error(calc_perm_p_value(statistic = Inf, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'statistic' failed: Must be finite.")

  # 'distribution'
  testthat::expect_error(calc_perm_p_value(statistic = statistic, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "'distribution' is missing.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = as.character(distribution), m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'distribution' failed: Must be of type 'numeric', not 'character'.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = c(NA, distribution), m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'distribution' failed: Contains missing values.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = NULL, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'distribution' failed: Must be of type 'numeric', not 'NULL'.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = c(distribution, Inf), m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'distribution' failed: Must be finite.")

  # 'm'
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "'m' is missing.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = "1", n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'm' failed: Must be of type 'count', not 'character'.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = c(1, 2), n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'm' failed: Must have length 1.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = NA, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'm' failed: May not be NA.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = -1, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'm' failed: Must be >= 1.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = NULL, n = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'm' failed: Must be of type 'count', not 'NULL'.")

  # 'n'
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "'n' is missing.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = "1", randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'n' failed: Must be of type 'count', not 'character'.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = c(1, 2), randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'n' failed: Must have length 1.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = NA, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'n' failed: May not be NA.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = -1, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'n' failed: Must be >= 1.")
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = NULL, randomization = FALSE, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'n' failed: Must be of type 'count', not 'NULL'.")

  # 'randomization'
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, n.rep = 10000, alternative = "two.sided"),
                         regexp = "'randomization' is missing.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = NA, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'randomization' failed: May not be NA.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = 1, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'randomization' failed: Must be of type 'logical flag', not 'double'.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = c(TRUE, TRUE), n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'randomization' failed: Must have length 1.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, NULL, n.rep = 10000, alternative = "two.sided"),
                         regexp = "Assertion on 'randomization' failed: Must be of type 'logical flag', not 'NULL'.",
                         fixed = TRUE)

  # 'n.rep'
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, alternative = "two.sided"),
                         regexp = "'n.rep' is missing.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = "10000", alternative = "two.sided"),
                         regexp = "Assertion on 'n.rep' failed: Must be of type 'count', not 'character'.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = c(1, 1000), alternative = "two.sided"),
                         regexp = "Assertion on 'n.rep' failed: Must have length 1.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = NULL, alternative = "two.sided"),
                         regexp = "Assertion on 'n.rep' failed: Must be of type 'count', not 'NULL'.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = -1, alternative = "two.sided"),
                         regexp = "Assertion on 'n.rep' failed: Must be >= 1.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = Inf, alternative = "two.sided"),
                         regexp = "Assertion on 'n.rep' failed: Must be of type 'count', not 'double'.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = TRUE, n.rep = 1000, alternative = "two.sided"),
                         regexp = "'n.rep' may not be larger than 252, the number of all splits.",
                         fixed = TRUE)

  # 'alternative'
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000),
                         regexp = "'alternative' is missing.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = 1),
                         regexp = "Assertion on 'alternative' failed: Must be element of set {'two.sided','greater','less'}, but types do not match (numeric != character).",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = c("greater", "two.sided")),
                         regexp = "Assertion on 'alternative' failed: Must be element of set {'two.sided','greater','less'}, but is not atomic scalar.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = "larger"),
                         regexp = "Assertion on 'alternative' failed: Must be element of set {'two.sided','greater','less'}, but is 'larger'.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = NA),
                         regexp = "Assertion on 'alternative' failed: Must be element of set {'two.sided','greater','less'}, but is 'NA'.",
                         fixed = TRUE)
  testthat::expect_error(calc_perm_p_value(statistic = statistic, distribution = distribution, m = 5, n = 5, randomization = FALSE, n.rep = 10000, alternative = NULL),
                         regexp = "Assertion on 'alternative' failed: Must be a subset of {'two.sided','greater','less'}, not 'NULL'.",
                         fixed = TRUE)

  ## Check output ----

  ## The output should be a numeric scalar

  # Permutation distribution
  checkmate::expect_number(calc_perm_p_value(statistic = statistic, distribution = distribution, m = m, n = n, randomization = FALSE, n.rep = 10000, alternative = "two.sided"), lower = 0, upper = 1)
  checkmate::expect_number(calc_perm_p_value(statistic = statistic, distribution = distribution, m = m, n = n, randomization = FALSE, n.rep = 10000, alternative = "greater"), lower = 0, upper = 1)
  checkmate::expect_number(calc_perm_p_value(statistic = statistic, distribution = distribution, m = m, n = n, randomization = FALSE, n.rep = 10000, alternative = "less"), lower = 0, upper = 1)
  testthat::expect_equal(calc_perm_p_value(statistic = statistic, distribution = distribution, m = m, n = n, randomization = FALSE, n.rep = 10000, alternative = "two.sided"), expected = 203/252)
  testthat::expect_equal(calc_perm_p_value(statistic = statistic, distribution = distribution, m = m, n = n, randomization = FALSE, n.rep = 10000, alternative = "greater"), expected = 203/252)
  testthat::expect_equal(calc_perm_p_value(statistic = statistic, distribution = distribution, m = m, n = n, randomization = FALSE, n.rep = 10000, alternative = "less"), expected = 50/252)

  # Randomization distribution
  checkmate::expect_number(calc_perm_p_value(statistic = statistic, distribution = distribution, m = m, n = n, randomization = TRUE, n.rep = 250, alternative = "two.sided"), lower = 0, upper = 1)
  checkmate::expect_number(calc_perm_p_value(statistic = statistic, distribution = distribution, m = m, n = n, randomization = TRUE, n.rep = 250, alternative = "greater"), lower = 0, upper = 1)
  checkmate::expect_number(calc_perm_p_value(statistic = statistic, distribution = distribution, m = m, n = n, randomization = TRUE, n.rep = 250, alternative = "less"), lower = 0, upper = 1)
})
