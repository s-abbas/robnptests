## ----------------------------------------------------------------------------
## Add random noise if too many values in the samples are equal
## ----------------------------------------------------------------------------

#' @title Add random noise to remove ties
#'
#' @description
#' \code{wobble} adds noise from a continuous uniform distribution to the
#' observations to remove ties.
#'
#' @template x
#' @template y
#' @template check_wobble
#'
#' @details
#' If \code{check = TRUE} the function checks whether all values in the two numeric
#' input vectors are distinct. If so, it returns the original values, otherwise
#' the ties are removed by adding noise from a continuous uniform distribution
#' to all observations. If \code{check = FALSE}, it simply determines the number
#' of digits and adds uniform noise.
#'
#' More precisely, we determine the minimum number of digits \code{d_min} in the sample
#' and then add random numbers from the U[-0.5 10^(-\code{d_min}), 0.5 10^(-\code{d_min})]
#' distribution to each of the observations.
#'
#' @return
#' A named list of length two containing the modified input samples \code{x} and
#' \code{y}.
#'
#' @references
#' \insertRef{FriGat07rank}{robnptests}
#'
#' @examples
#' x <- rnorm(20); y <- rnorm(20); x <- round(x)
#' wobble(x, y)
#'
#' @export

wobble <- function(x, y, check = TRUE) {

  # Check input arguments ----
  stopifnot("'x' is missing" = !missing(x),
            "'y' is missing" = !missing(y)
  )

  checkmate::assert_double(x, any.missing = FALSE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_double(y, any.missing = FALSE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_flag(check, na.ok = FALSE, null.ok = FALSE)

  # Wobbling ----
  # Determine number of different values in both samples and in joint sample
  no.values.x <- length(unique(x))
  no.values.y <- length(unique(y))
  no.values.xy <- length(unique(c(x, y)))

  if (check) {
    # Determine number of different values in both samples and in joint sample
    no.values.x <- length(unique(x))
    no.values.y <- length(unique(y))
    no.values.xy <- length(unique(c(x, y)))
    if (no.values.x == length(x) & no.values.y == length(y) & no.values.xy == length(c(x, y))) {
    # If all values are distinct, return original observations
     return(list(x = x, y = y))
    }
  }

  z <- c(x, y)

  # Maximal number of digits after decimal point:
  # this expression always returns two values (part before and part after the decimal
  # point)
  # max.digits <- max(sapply(strsplit(as.character(z), "\\."), length)) - 1
  # Instead:
  digits <- as.vector(sapply(as.character(z),
                             function(x) nchar(unlist(strsplit(x, "\\."))[2])))
  # if the values are discrete, 'strsplit' returns NA
  digits[is.na(digits)] <- 0

  z.wobble <- z + stats::runif(length(z),
                        -0.5*10^(-min(digits)), 0.5*10^(-min(digits)))

  return(list(x = z.wobble[1:length(x)], y = z.wobble[(length(x) + 1):length(c(x, y))]))
}
