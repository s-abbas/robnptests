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
  ## Error handling
  if (gamma < 0 || gamma > 0.5) {
    stop ("gamma has to be in [0, 0.5]")
  }

  ## NA handling
  if (!na.rm & any(is.na(x))) {
    return(NA)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Calculate winsorized variance
  n <- length(x)

  r <- floor(gamma * n)

  x.sort <- sort(x)
  x.lower <- x.sort[r + 1]
  x.upper <- x.sort[n - r]
  x.sort[x.sort < x.lower] <- x.lower
  x.sort[x.sort > x.upper] <- x.upper

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
#'
#' @details
#' For definitions of the scale estimators see Fried and Dehling (2011).
#'
#'
#' @return
#' An estimate of the pooled variance of the two samples.
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{FriDeh11robu}{robTests}
#'
#'
#' @export

rob_var <- function(x, y, na.rm = FALSE, type = c("S1", "S2", "S3", "S4")) {
  type <- match.arg(type)

  ## Error handling
  if (!(type %in% c("S1", "S2", "S3", "S4"))) {
    stop("'type' needs to be one of 'S1', 'S2', 'S3', 'S4'.")
  }
  stopifnot(
    is.numeric(x),
    is.numeric(y)
  )

  if ((length(unique(x)) == 1) | (length(unique(y)) == 1)) {
    warning("At least one of the input vectors is constant.", call. = FALSE)
  }

  ## NA handling
  if (!na.rm & (any(is.na(x)) || any(is.na(y)))) {
    return(NA)
  } else if (na.rm & (any(is.na(x)) || any(is.na(y)))) {
    x <- as.vector(stats::na.omit(x))
    y <- as.vector(stats::na.omit(y))
  }

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
    stop("Estimate of scale is 0 although the data is not constant. Consider using a different estimator or setting wobble = TRUE in the function call.",
            call. = FALSE)
  }

  return(est)
}
