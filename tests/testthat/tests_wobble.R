context("wobbling")

testthat::test_that("wobbling works correctly", {

  x <- c(1, 2, 3, 4, 4)
  y <- c(6, 7, 8, 9, 9)

  testthat::expect_error(wobble(x = x), regexp = "'y' is missing", fixed = TRUE)
  testthat::expect_error(wobble(y = y), regexp = "'x' is missing", fixed = TRUE)

  # Check what happens if a character is handed over
  testthat::expect_error(wobble(x = x, y = "character"),
                         regexp = "Assertion on 'y' failed: Must be of type 'double', not 'character'.")
  testthat::expect_error(wobble(x = "character", y = y),
                         regexp = "Assertion on 'x' failed: Must be of type 'double', not 'character'.")
  testthat::expect_error(wobble(x = "character", y = "character"),
                         regexp = "Assertion on 'x' failed: Must be of type 'double', not 'character'.")

  # Check for NAs
  testthat::expect_error(wobble(x = c(x, NA), y = y),
                         regexp = "Assertion on 'x' failed: Contains missing values (element 6).",
                         fixed = TRUE)
  testthat::expect_error(wobble(x = x, y = c(y, NA)),
                         regexp = "Assertion on 'y' failed: Contains missing values (element 6).",
                         fixed = TRUE)

  # Check that check works
  testthat::expect_error(wobble(x = x, y = y, check = "character"),
                         regexp = "Assertion on 'check' failed: Must be of type 'logical flag', not 'character'.")
  testthat::expect_error(wobble(x = x, y = y, check = NA),
                         regexp = "Assertion on 'check' failed: May not be NA.")
  testthat::expect_error(wobble(x = x, y = y, check = 3),
                         regexp = "Assertion on 'check' failed: Must be of type 'logical flag', not 'double'.")
  testthat::expect_error(wobble(x = x, y = y, check = c(TRUE, NA, 3)),
                         regexp = "Assertion on 'check' failed: Must be of type 'logical flag', not 'double'.")
  testthat::expect_error(wobble(x = x, y = y, check = c(TRUE, FALSE)),
                         regexp = "Assertion on 'check' failed: Must have length 1.")

  # Now on the actual function:
  # The function should only wobble if there are bindings in the sample
  testthat::expect_equal(wobble(x = unique(x), y = unique(y), check = TRUE)$x, unique(x))
  testthat::expect_equal(wobble(x = unique(x), y = unique(y), check = TRUE)$y, unique(y))

  # Check whether the wobbled samples are reproducible if a seed is set
  testthat::expect_equal({set.seed(1212); wobble(x = x, y = y, check = TRUE)$x},
                         {set.seed(1212); wobble(x = x, y = y, check = TRUE)$x})
  testthat::expect_equal({set.seed(1212); wobble(x = x, y = y, check = TRUE)$y},
                         {set.seed(1212); wobble(x = x, y = y, check = TRUE)$y})

  # Check whether the function also wobbles if no bindings are present, but check = FALSE
  # i.e. the samples should change!
  testthat::expect_false(isTRUE(all.equal(wobble(x = x, y = y, check = TRUE)$y, y)))
  testthat::expect_false(isTRUE(all.equal(wobble(x = x, y = y, check = TRUE)$x, x)))

  # Check whether the wobbling really takes place in d_min
  # Help function:  (IF ANY CHANGES ARE MADE IN wobble(), THIS MIGHT NEED ADJUSTING)
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
  testthat::expect_equal(get_dmin(x+ 0.1, y + 0.1), 1)

}
)
