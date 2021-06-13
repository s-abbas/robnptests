testthat::test_that("wobbling works correctly", {

  # Exemplary input vectors ----
  x <- c(1, 2, 3, 4, 4)
  y <- c(6, 7, 8, 9, 9)

  # Wobbling only if there are bindings in the sample ----
  testthat::expect_equal(wobble(x = unique(x), y = unique(y), check = TRUE)$x, unique(x))
  testthat::expect_equal(wobble(x = unique(x), y = unique(y), check = TRUE)$y, unique(y))

  # Samples should be reproducible if seed is set ----
  testthat::expect_equal({set.seed(1212); wobble(x = x, y = y, check = TRUE)$x},
                         {set.seed(1212); wobble(x = x, y = y, check = TRUE)$x})
  testthat::expect_equal({set.seed(1212); wobble(x = x, y = y, check = TRUE)$y},
                         {set.seed(1212); wobble(x = x, y = y, check = TRUE)$y})

  # If 'check' = FALSE, ----
  # wobbling should be performed even if no bindings are present,
  # i.e. the samples should change!
  set.seed(1212)
  testthat::expect_false(isTRUE(all.equal(wobble(x = x, y = y, check = TRUE)$y, y)))
  set.seed(1212)
  testthat::expect_false(isTRUE(all.equal(wobble(x = x, y = y, check = TRUE)$x, x)))

  # Wobbling should take place in decimal place 'd_min' ----
  # Helper function: (if changes are made in 'wobble', this might need adjustment)
  get_dmin <- function(x, y) {
    z <- c(x, y)
    digits <- as.vector(sapply(as.character(z),
                               function(x) nchar(unlist(strsplit(x, "\\."))[2])))
    digits[is.na(digits)] <- 0
    return(min(digits))
  }

  testthat::expect_equal(get_dmin(x, y), 0)
  testthat::expect_equal(get_dmin(c(x, 1.1), y), 0)
  testthat::expect_equal(get_dmin(x, c(y, 1.2, 1.23)), 0)
  testthat::expect_equal(get_dmin(x + 0.1, y + 0.1), 1)

}
)
