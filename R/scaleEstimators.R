## ----------------------------------------------------------------------------
## Different scale estimators for one and two samples
## ----------------------------------------------------------------------------

#' @title Winsorized variance
#'
#' @description
#' \code{win_var} calculates the winsorized variance of a sample.
#'
#'
#' @template x
#' @template gamma_winsorized
#' @template na_rm
#'
#' @return A list containing the following items:
#' \item{var}{winsorized variance.}
#' \item{h}{degrees of freedom used for tests based on trimmed means and the
#' winsorized variance.}
#'
#'
#' @examples
#' ## Generate random sample
#' set.seed(108)
#' x <- rnorm(10)
#'
#' ## Compute 20% winsorized variance
#' win_var(x, gamma = 0.2)
#'
#' @export

win_var <- function(x, gamma = 0, na.rm = FALSE) {

  ## Check input arguments
  stopifnot("'x' is missing." = !missing(x))

  checkmate::assert_numeric(x, min.len = 2, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_number(gamma, lower = 0, upper = 0.5, na.ok = FALSE, finite = TRUE, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)

  ## Remove missing values in 'x' ----
  if (!na.rm & any(is.na(x))) {
    return(list(var = NA_real_, h = NA_real_))
  } else if (na.rm) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Error message if all values in 'x' are equal ----
  if (length(unique(x)) == 1) {
    stop("All values in '", deparse(substitute(x)), "' are equal. The scale estimate is '0' and the test statistic cannot be computed.")
  }

  ## Calculate winsorized variance ----
  n <- length(x)

  r <- floor(gamma * n)

  x.sort <- sort(x)
  x.lower <- x.sort[r + 1]
  x.upper <- x.sort[n - r]
  x.sort[which(x.sort < x.lower)] <- x.lower
  x.sort[which(x.sort > x.upper)] <- x.upper

  res <- 1 / (n - 1) * sum((x.sort - mean(x.sort)) ^ 2)
  h <- n - 2 * r

  return(list(var = res, h = h))
}

#' @title Robust scale estimators based on median absolute deviation
#'
#' @description
#' \code{rob_var} calculates a variance estimator for the within-sample variance
#' based on two samples.
#'
#' @template x
#' @template y
#' @template scale_type
#' @template na_rm
#' @template check_for_zero
#'
#' @details
#' For definitions of the scale estimators see Fried and Dehling (2011).
#'
#' If \code{check.for.zero = TRUE}, a warning is printed when the scale estimate
#' is zero. This argument is only included as the function is used in
#' \code{\link{rob_perm_statistic}} to computed values of robust test statistics
#' where the scale estimate is a standardization. A scale estimate of zero leads
#' to a non-existing test statistic, so that the corresponding test cannot be
#' performed.
#'
#' @return
#' An estimate of the pooled variance of the two samples.
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{FriDeh11robu}{robnptests}
#'
#' @export

rob_var <- function(x, y, type = c("S1", "S2", "S3", "S4"), na.rm = FALSE, check.for.zero = FALSE) {

  ## Check input arguments
  stopifnot("'x' is missing." = !missing(x))
  stopifnot("'y' is missing." = !missing(y))

  checkmate::assert_numeric(x, min.len = 2, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(y, min.len = 2, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_subset(type, choices = c("S1", "S2", "S3", "S4"), empty.ok = FALSE)
  checkmate::assert_character(type, min.chars = 1, ignore.case = FALSE, all.missing = FALSE, min.len = 1, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)
  checkmate::assert_flag(check.for.zero, na.ok = FALSE, null.ok = FALSE)

  ## Match 'type' ----
  type <- match.arg(type)

  # Remove missing values in 'x' and 'y' ----
  if (!na.rm & (any(is.na(x)) || any(is.na(y)))) {
    return(NA_real_)
  } else if (na.rm) {
    x <- as.vector(stats::na.omit(x))
    y <- as.vector(stats::na.omit(y))
  }

  ## Compute scale estimates ----
  if (type == "S1") {
    xcomb <- utils::combn(x, 2)
    ycomb <- utils::combn(y, 2)
    xabs <- abs(xcomb[1, ] - xcomb[2, ])
    yabs <- abs(ycomb[1, ] - ycomb[2, ])
    est <- stats::median(c(xabs, yabs))
  } else if (type == "S2") {
    z <- c(x - stats::median(x), y - stats::median(y))
    zcomb <- utils::combn(z, 2)
    est <- stats::median(abs(zcomb[1, ] - zcomb[2, ]))
  } else  if (type == "S3") {
    est <- 2 * stats::median(c(abs(x - stats::median(x)), abs(y - stats::median(y))))
  } else if (type == "S4") {
    est <- stats::median(abs(x - stats::median(x)) + stats::median(abs(y - stats::median(y))))
  }

  ## Error message, if 'est' = 0 ----
  if (check.for.zero & est == 0) {
    if (length(unique(c(x, y))) == 1) {
      # All values in 'x' and 'y' are equal
      stop("All values in '", deparse(substitute(x)), "' ", "and '", deparse(substitute(y)), "' ", "are equal. The scale estimate is zero and the test statistic cannot be computed.", call. = FALSE)
    } else {
      # Scale estimate is zero although data are not constant
      stop("A scale estimate of zero occured although the data are not constant. Consider using a different scale estimator or set 'wobble = TRUE' in the function call.", call. = FALSE)
    }
  }

  return(est)
}
