# Trimmed t ----
testthat::test_that("trimmed_t works correctly", {

  # Exemplary input vectors ----
  x <- c(22, 23, 31, 71, 74)
  y <- c(4, 4, 8, 19, 26, 74)

  # Results for 'x' computed by hand with 'gamma = 0.2':
  # trimmed_mean(x) = 41.66667
  # win_var(x) = 627.2
  # h = 3
  #
  # Results for 'y' computed by hand with 'gamma = 0.2':
  # trimmed_mean(y) = 14.25
  # win_var(y) = 109.5
  # h = 4
  #
  # Test statistic: 1.451922

  # Removal of missing values ----
  testthat::expect_equal(trimmed_t(x = c(x, NA), y = y, na.rm = TRUE),
                         trimmed_t(x = x, y = y))

  testthat::expect_equal(trimmed_t(x = c(x, NA), y = y, na.rm = FALSE),
                         list(statistic = NA_real_,
                              estimates = c(NA_real_, trim_mean(y)),
                              df = NA_real_))

  # Location and scale invariance of the test statistic ----
  testthat::expect_equal(trimmed_t(x = 2 * x + 3, y = 2 * y + 3)$statistic,
                         trimmed_t(x = x, y = y)$statistic)

  # Output check ----

  # The output should be a list of numeric scalars
  checkmate::expect_list(trimmed_t(x = x, y = y), types = rep("numeric", 3))

  checkmate::expect_list(trimmed_t(x = c(NA, x), y = c(NA, y)),
                         types = rep("numeric", 3))

  # The output should be equal to the manually computed values
  testthat::expect_equal(trimmed_t(x = x, y = y, gamma = 0.2),
                         list(statistic = 1.451922,
                              estimates = c(41.66667, 14.25),
                              df = 5),
                         tolerance = 1e-07)

})

# Test statistics of robust permutation tests ----
testthat::test_that("rob_perm_statistic works correctly", {

  # Exemplary input vectors ----
  x <- c(22, 23, 31, 71, 74)
  y <- c(4, 4, 8, 19, 26, 74)

  # Results for 'x' computed by hand:
  # hodges_lehmann(x) =  47.5
  # median(x) = 31
  #
  # Results for 'y' computed by hand:
  # hodges_lehmann(y) = 15
  # median(y) = 13.5
  #
  # Results for HL2-estimator computed by hand:
  # hodges_lehmann_2sample(x, y) = 18.5
  #
  # Variance for 'S1':
  # rob_var(x, y, type = "S1") = 22
  #
  # Variance for 'S2':
  # rob_var(x, y, type = "S2") = 20.5
  #
  # Variance for 'S3':
  # rob_var(x, y, type = "S3") = 19
  #
  # Variance for 'S4':
  # rob_var(x, y, type = "S3") = 18.5
  #
  # Test statistic for 'type = HL11': 1.477273
  # Test statistic for 'type = HL12': 1.585366
  # Test statistic for 'type = HL21': 0.8409091
  # Test statistic for 'type = HL22': 0.902439
  # Test statistic for 'type = MED1': 0.9210526
  # Test statistic for 'type = MED2': 0.9459459

  # Removal of missing values ----
  testthat::expect_equal(rob_perm_statistic(x = c(x, NA), y = y, na.rm = TRUE),
                         rob_perm_statistic(x = x, y = y))

  testthat::expect_equal(rob_perm_statistic(x = c(x, NA), y = y, na.rm = FALSE,
                                            type = "HL11"),
                         list(statistic = NA_real_,
                              estimates = c(NA_real_, hodges_lehmann(y))))

  testthat::expect_equal(rob_perm_statistic(x = c(x, NA), y = y, na.rm = FALSE,
                                            type = "HL21"),
                         list(statistic = NA_real_,
                              estimates = c(NULL, NULL)))

  # Location and scale invariance of the test statistic ----
  types <- c("HL11", "HL12", "HL21", "HL22", "MED1", "MED2")

  for (i in seq_along(types)) {
    testthat::expect_equal(rob_perm_statistic(x = 2 * x + 3, y = 2 * y + 3,
                                              type = types[i])$statistic,
                           rob_perm_statistic(x = x, y = y,
                                              type = types[i])$statistic)
  }

  # Output check ----

  # The output should be a list of numeric scalars
  for (i in seq_along(types)) {
    if (!(types[i] %in% c("HL21", "HL22"))) {
      checkmate::expect_list(rob_perm_statistic(x = x, y = y, type = types[i]),
                             types = rep("numeric", 3))
    } else {
      checkmate::expect_list(rob_perm_statistic(x = x, y = y, type = types[i]),
                             types = c("numeric", "numeric", "NULL"))
    }

    if (!(types[i] %in% c("HL21", "HL22"))) {
      checkmate::expect_list(rob_perm_statistic(x = c(NA, x), y = c(NA, y),
                                                type = types[i]),
                             types = rep("numeric", 3))
    } else {
      checkmate::expect_list(rob_perm_statistic(x = c(NA, x), y = c(NA, y),
                                                type = types[i]),
                             types = c("numeric", "numeric", "NULL"))
    }
  }

  # The output should be equal to the manually computed values

  # HL1-statistics
  testthat::expect_equal(rob_perm_statistic(x = x, y = y, type = "HL11"),
                         list(statistic = 1.477273,
                              estimates = c(47.5, 15)),
                         tolerance = 1e-06)

  testthat::expect_equal(rob_perm_statistic(x = x, y = y, type = "HL12"),
                         list(statistic = 1.585366,
                              estimates = c(47.5, 15)),
                         tolerance = 1e-06)

  # HL2-statistics
  testthat::expect_equal(rob_perm_statistic(x = x, y = y, type = "HL21"),
                         list(statistic = 0.8409091,
                              estimates = NULL),
                         tolerance = 1e-06)

  testthat::expect_equal(rob_perm_statistic(x = x, y = y, type = "HL22"),
                         list(statistic = 0.902439,
                              estimates = NULL),
                         tolerance = 1e-06)

  # MED-statistics
  testthat::expect_equal(rob_perm_statistic(x = x, y = y, type = "MED1"),
                         list(statistic = 0.9210526,
                              estimates = c(31, 13.5)),
                         tolerance = 1e-06)

  testthat::expect_equal(rob_perm_statistic(x = x, y = y, type = "MED2"),
                         list(statistic = 0.9459459,
                              estimates = c(31, 13.5)),
                         tolerance = 1e-06)
})

# M-tests ----
testthat::test_that("m_test_statistic works correctly", {
  psi.funs <- c("huber", "hampel", "bisquare")

  # Exemplary input vectors ----
  x <- c(22, 23, 31, 71, 74)
  y <- c(4, 4, 8, 19, 26, 74)

  # Removal of missing values ----
  for (i in seq_along(psi.funs)) {
    suppressWarnings(
      testthat::expect_equal(m_test_statistic(x = c(x, NA), y = y,
                                              psi = psi.funs[i]),
                             m_test_statistic(x = x, y = y, psi = psi.funs[i]))
    )

    suppressWarnings(
      testthat::expect_error(m_test_statistic(x = c(x[1:2], NA, NA, NA), y = y,
                                              psi = psi.funs[i]))
    )

    # Location and scale invariance of the test statistic ----
    testthat::expect_equal(m_test_statistic(x = 2 * x + 3, y = 2 * y + 3,
                                            psi = psi.funs[i])$statistic,
                           m_test_statistic(x = x, y = y,
                                            psi = psi.funs[i])$statistic)

    # Output check ----

    # The output should be a list of numeric scalars
    checkmate::expect_list(trimmed_t(x = x, y = y), types = rep("numeric", 2))

    testthat::expect_warning(m_test_statistic(x = c(x, NA), y = y,
                                              psi = psi.funs[i]))

    testthat::expect_warning(m_test_statistic(x = x, y = c(y, NA),
                                              psi = psi.funs[i]))

    # Create and compare snapshots of test output
    # The output of the function cannot be computed manually
    testthat::expect_snapshot_output(m_test_statistic(x = x, y = y,
                                                      psi = psi.funs[i]))
  }
})
