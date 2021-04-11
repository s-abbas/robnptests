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
#' @template gamma_winsorized_variance
#' @template na_rm
#'
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
#'
#' @details
#' For definitions of the scale estimators see Fried and Dehling (2011).
#'
#' @return
#' An estimate of the pooled variance of the two samples.
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{FriDeh11robu}{robTests}
#'
#' @export

rob_var <- function(x, y, type = c("S1", "S2", "S3", "S4"), na.rm = FALSE) {

  ## Check input arguments
  stopifnot("'x' is missing." = !missing(x))
  stopifnot("'y' is missing." = !missing(y))

  checkmate::assert_numeric(x, min.len = 2, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(y, min.len = 2, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_subset(type, choices = c("S1", "S2", "S3", "S4"), empty.ok = FALSE)
  checkmate::assert_character(type, min.chars = 1, ignore.case = FALSE, all.missing = FALSE, min.len = 1, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)

  ## Match 'type' ----
  type <- match.arg(type)

  # Remove missing values in 'x' and 'y' ----
  if (!na.rm & (any(is.na(x)) || any(is.na(y)))) {
    return(NA_real_)
  } else if (na.rm) {
    x <- as.vector(stats::na.omit(x))
    y <- as.vector(stats::na.omit(y))
  }

  ## Error message if all values in 'x' are equal ----
  if (length(unique(x)) == 1 & length(unique(y)) == 1) {
    stop("All values in '", deparse(substitute(x)), "' ", "and '", deparse(substitue(y)), "' ", "are equal. The scale estimate is '0' and the test statistic cannot be computed.")
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

  if (est == 0 & (length(unique(x)) > 1 & length(unique(y)) > 1)) {
    stop( "Estimate of scale is '0' although the data are not constant.
          Consider using a different estimator or setting wobble = TRUE in the function call. Otherwise, the test statistic cannot be computed.")
  }

  return(est)
}
