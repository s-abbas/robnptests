## ----------------------------------------------------------------------------
## Add random noise if too many values in the samples are equal
## ----------------------------------------------------------------------------

#' @title Function that undiscretizes rounded samples
#'
#' @description
#' \code{wobble} makes a discrete sample with duplicated values continuous by adding uniform noise to the
#' discrete observations
#'
#' @template x
#' @template y
#' @param check logical value indicating whether the samples should be checked for bindings prior to
#'              adding uniform noise or not, defaults to \code{TRUE}
#'
#' @details
#' If \code{check = TRUE} the function checks whether all values in the two numeric input vectors are distinct.
#' If so, it returns the original values, otherwise the values are made continuous by adding uniform
#' noise. If \code{check = FALSE}, it simply determines the number of digits and adds uniform noise in order to
#' make the sample "more" continuous
#'
#' Precisely, we determine the minimum number of digits d_min in the sample
#' and then add random variables from the U[-0.5 10^(-d_min), 0.5 10^(-d_min)] distribution to each
#' of the observations.
#'
#' @return
#' A list of length two containing the modified \code{x} and \code{y}.
#'
#' @references
#' \insertRef{FriGat07rank}{robTests}
#'
#' @examples
#' x <- rnorm(20); y <- rnorm(20); x <- round(x)
#' wobble(x, y)
#'
#'
#' @export


wobble <- function(x, y, check = TRUE) {

  # Check input arguments:
  checkmate::assert_double(x, any.missing = FALSE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_double(y, any.missing = FALSE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_flag(check, na.ok = FALSE, null.ok = FALSE)


  ## Determine number of different values in both samples and in joint sample
  no.values.x <- length(unique(x))
  no.values.y <- length(unique(y))
  no.values.xy <- length(unique(c(x, y)))

  if (check) {
    ## Determine number of different values in both samples and in joint sample
    no.values.x <- length(unique(x))
    no.values.y <- length(unique(y))
    no.values.xy <- length(unique(c(x, y)))
    if (no.values.x == length(x) & no.values.y == length(y) & no.values.xy == length(c(x, y))) {
    ## If all values are distinct, return original observations
     return(list(x = x, y = y))
    }
  }

  z <- c(x, y)

  ## Maximal number of digits after decimal point:
  # this expression always returns two (part before and part after the decimal
  # point)
  # max.digits <- max(sapply(strsplit(as.character(z), "\\."), length)) - 1
  # Instead:
  digits <- as.vector(sapply(as.character(z),
                             function(x) nchar(unlist(strsplit(x, "\\."))[2])))
  digits[is.na(digits)] <- 0 # if the values are discrete we will get NA from the strsplit

  # I propose we use the full range of values between the observations:
  ## Add random noise from a uniform distribution
  z.wobble <- z + stats::runif(length(z),
                        -0.5*10^(-min(digits)), 0.5*10^(-min(digits)))

  #z.wobble <- z + runif(length(z), min = -10^(-(max.digits + 1)), max = 10^(-(max.digits + 1)))

  return(list(x = z.wobble[1:length(x)], y = z.wobble[(length(x) + 1):length(c(x, y))]))
}





