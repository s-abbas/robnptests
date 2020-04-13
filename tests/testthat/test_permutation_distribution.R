context("Permutation distribution")

testthat::test_that("perm_distribution works correctly", {

  ## ___________________________________________________________________________
  ## Number of randomly drawn splits for randomization tests must not be larger
  ## than the number of all splits
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(5)
  y <- rnorm(5)

  testthat::expect_error(perm_distribution(x, y, type = "D1S1", sampled = TRUE, n.rep = 1000))

  ## ___________________________________________________________________________
  ## Check length of vector with permutation distribution
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(5)
  y <- rnorm(5)

  distribution <- perm_distribution(x, y, type = "D1S1", sampled = FALSE)

  testthat::expect_equal(length(distribution), choose(length(x) + length(y), length(x)))

  ## ___________________________________________________________________________
  ## Check length of vector with randomization distribution
  ## ___________________________________________________________________________
  set.seed(108)

  x <- rnorm(10)
  y <- rnorm(10)

  distribution <- perm_distribution(x, y, type = "D1S1", sampled = TRUE, n.rep = 1000)

  testthat::expect_equal(length(distribution), 1000)
})
