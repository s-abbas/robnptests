
## Computation of permutation distribution ----
testthat::test_that("perm_distribution works correctly", {

  # testthat::skip_on_cran()

  ## Generate exemplary input vectors
  x <- 1:5
  y <- 6:10
  m <- length(x)
  n <- length(y)

  ## Checks for input arguments ----
  # 'x'
  testthat::expect_error(perm_distribution(y = y, type = "HL11"),
                         regexp = "'x' is missing.")
  testthat::expect_error(perm_distribution(x = c(x, NA), y = y, type = "HL11"),
                         regexp = "Assertion on 'x' failed: Contains missing values.")
  testthat::expect_error(perm_distribution(x = NULL, y = y, type = "HL11"),
                         regexp = "Assertion on 'x' failed: Must be of type 'numeric', not 'NULL'.")
  testthat::expect_error(perm_distribution(x = c(x, Inf), y = y, type = "HL11"),
                         regexp = "Assertion on 'x' failed: Must be finite.")
  testthat::expect_error(perm_distribution(x = x[1:2], y = y, type = "HL11"),
                         regexp = "Assertion on 'x' failed: Must have length >= 5, but has length 2.")

  # 'y'
  testthat::expect_error(perm_distribution(x = x, type = "HL22"),
                         regexp = "'y' is missing.")
  testthat::expect_error(perm_distribution(x = x, y = c(NA, y), type = "HL22"),
                         regexp = "Assertion on 'y' failed: Contains missing values.")
  testthat::expect_error(perm_distribution(x = x, y = NULL, type = "HL11"),
                         regexp = "Assertion on 'y' failed: Must be of type 'numeric', not 'NULL'.")
  testthat::expect_error(perm_distribution(x = x, y = c(-Inf, y), type = "HL11"),
                         regexp = "Assertion on 'y' failed: Must be finite.")
  testthat::expect_error(perm_distribution(x = x, y = y[3:5], type = "HL11"),
                         regexp = "Assertion on 'y' failed: Must have length >= 5, but has length 3.")

  # 'type'
  testthat::expect_error(perm_distribution(x = x, y = y),
                         regexp = "'type' is missing.")
  testthat::expect_error(perm_distribution(x = x, y = y, type = 1),
                         regexp = "Assertion on 'type' failed: Must be element of set {'HL11','HL12','HL21','HL22','MED1','MED2'}, but types do not match (numeric != character).", fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = c("HL11", "HL12")),
                         regexp = "Assertion on 'type' failed: Must be element of set {'HL11','HL12','HL21','HL22','MED1','MED2'}, but is not atomic scalar.", fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = "HL23"), regexp = "Assertion on 'type' failed: Must be element of set {'HL11','HL12','HL21','HL22','MED1','MED2'}, but is 'HL23'.", fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = NA), regexp = "Assertion on 'type' failed: Must be element of set {'HL11','HL12','HL21','HL22','MED1','MED2'}, but is 'NA'.", fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = NULL), regexp = "Assertion on 'type' failed: Must be a subset of {'HL11','HL12','HL21','HL22','MED1','MED2'}, not 'NULL'.", fixed = TRUE)

  # 'randomization'
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = NA),
                         regexp = "Assertion on 'randomization' failed: May not be NA.",
                         fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = 1),
                         regexp = "Assertion on 'randomization' failed: Must be of type 'logical flag', not 'double'.",
                         fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = c(TRUE, TRUE)),
                         regexp = "Assertion on 'randomization' failed: Must have length 1.",
                         fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = NULL),
                         regexp = "Assertion on 'randomization' failed: Must be of type 'logical flag', not 'NULL'.",
                         fixed = TRUE)

  # 'n.rep'
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = FALSE, n.rep = NA),
                         regexp = "Assertion on 'n.rep' failed: May not be NA.",
                         fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = FALSE, n.rep = "10000"),
                         regexp = "Assertion on 'n.rep' failed: Must be of type 'count', not 'character'.",
                         fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = FALSE, n.rep = c(1, 100)),
                         regexp = "Assertion on 'n.rep' failed: Must have length 1.",
                         fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = FALSE, n.rep = NULL),
                         regexp = "Assertion on 'n.rep' failed: Must be of type 'count', not 'NULL'.",
                         fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = FALSE, n.rep = -100),
                         regexp = "Assertion on 'n.rep' failed: Must be >= 1.",
                         fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = FALSE, n.rep = Inf),
                         regexp = "Assertion on 'n.rep' failed: Must be of type 'count', not 'double'.",
                         fixed = TRUE)
  testthat::expect_error(perm_distribution(x = x, y = y, type = "MED1", randomization = TRUE, n.rep = 1000),
                         regexp = "'n.rep' may not be larger than 252, the number of all splits.",
                         fixed = TRUE)


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
  testthat::expect_error(m_est_perm_distribution(x = x, y = y, psi = "huber", k = 1.345),
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
