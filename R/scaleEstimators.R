## ----------------------------------------------------------------------------
## Different scale estimators for one and two samples
## ----------------------------------------------------------------------------

#' @title Winsorized variance
#'
#' @description
#' \code{win_var} calculates the winsorized variance of a sample.
#'
#' @inheritParams trim_mean
#'
#' @return A list containing the following items:
#' \item{var}{The winsorized variance, and}
#' \item{h}{the degrees of freedom used for tests based on trimmed means and the
#' winsorized variance.}
#'
#' @details
#'
#' @examples
#' x <- rnorm(10)
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
    x <- stats::na.omit(x)
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

#' @title Asymmetrically trimmed variance
#'
#' @description
#' \code{asym_win_var} calculates the asymmetrically trimmed variance using different
#' skewness selector statistics.
#'
#' @inheritParams trim_mean
#' @param type specifies the skewness selector statistic used for trimming/winsorizing,
#'  must be in \code{"Q2", "SK2"} and \code{"SK5"}. Default is \code{"Q2"}.
#'
#' @return A list containing the following items:
#' \item{var}{The asymmetrically trimmed variance, and}
#' \item{h}{the degrees of freedom used for tests based on  asymmetrically trimmed means and the
#' asymmetrically winsorized variance.}
#'
#' @examples
#' x <- rnorm(10)
#' asym_win_var(x, type = "SK5")
#'
#' @export

asym_win_var <- function(x, type = c("Q2", "SK2", "SK5"), na.rm = FALSE) {
  ## NA handling
  if (!na.rm & any(is.na(x))) {
    return(NA)
  } else if (na.rm & any(is.na(x))) {
    x <- stats::na.omit(x)
  }

  ## Sample size and ordered sample
  m <- length(x)
  x.sort <- sort(x)

  ## Calculate asymmetrically trimmed mean
  mean.x <- asym_trimmed_mean(x, type = type)

  if (type == "Q2") {
    if (floor(0.05 * m) == 0) {
      ## If sample size is too small for trimming, calculate ordinary
      ## empirical variance
      res <- stats::var(x)

      return(list(var = res, h = m))
    }

    gamma <- 0.1

    ## Determine lower trimming proportion
    U.005 <- mean(x.sort[(m - floor(0.05 * m) + 1):m])
    L.005 <- mean(x.sort[1:floor(0.05 * m)])
    U.05 <- mean(x.sort[(m - floor(0.5 * m) + 1):m])
    L.05 <- mean(x.sort[1:floor(0.5 * m)])

    gamma.lower <- gamma * (U.005 - L.005)/(U.005 - L.005 + U.05 - L.05)
  } else if (type == "SK2") {
    gamma <- 0.1

    gamma.lower <- gamma * (min(x) - stats::median(x))/(min(x) - max(x))
  } else if (type == "SK5") {
    gamma <- 0.25

    gamma.lower <- gamma * (min(x) - mean(x))/(min(x) - max(x))
  }

  ## Determine upper trimming proportion
  gamma.upper <- gamma - gamma.lower

  ## Replace observations in sorted sample to compute winsorized variance
  m.lower <- 1 + floor(gamma.lower * m + 0.5)
  m.upper <- m - floor(gamma.upper * m + 0.5)

  x.sort[1:(m.lower - 1)] <- x.sort[m.lower]
  x.sort[(m.upper + 1):m] <- x.sort[m.upper]

  ## Calculate variance
  res <- 1/(m.upper - m.lower + 1) * sum((x.sort - mean.x)^2)

  h <- m.upper - m.lower + 1

  return(list(var = res, h = h))
}


#' @title Robust 2 sample variance estimates based on medians
#'
#' @description
#' \code{rob_var} calculates a variance estimate based on median differences
#'
#' @param x numeric vector of observations.
#' @param y numeric vector of observations.
#' @param type in \code{"S1"}, \code{"S2"}, \code{"S3"} and \code{"S4"} for
#' respective variance estimator, see details for description of the scale estimators.
#' @param na.rm a logical value indicating whether NA values in x and y should be stripped before the computation proceeds.
#'
#' @details
#' For definitions of the scale estimators see Fried and Dehling (2008).
#'
#' @return
#' An estimate of the pooled variance of the two samples.
#'
#' @references
#' \insertRef{FriDeh11robu}{robTests}
#' @import utils
#' @export

rob_var <- function(x, y, na.rm = FALSE, type = c("S1", "S2", "S3", "S4")) {
  type <- match.arg(type)

  ## Error handling
  if (!(type %in% c("S1", "S2", "S3", "S4"))) {
    stop("type needs to be one of 'S1', 'S2', 'S3', 'S4'")
  }
  stopifnot(
    is.numeric(x),
    is.numeric(y)
  )

  ## NA handling
  if (!na.rm & any(is.na(x)) || any(is.na(y))) {
    return(NA)
  } else if (na.rm & any(is.na(x)) || is.na(y)) {
    x <- stats::na.omit(x)
    y <- stats::na.omit(y)
  }

  if (type == "S1") {
    xcomb <- utils::combn(x, 2)
    ycomb <- utils::combn(y, 2)
    xabs <- abs(xcomb[1, ] - xcomb[2, ])
    yabs <- abs(ycomb[1, ] - ycomb[2, ])
    return(stats::median(c(xabs, yabs)))

  } else if (type == "S2") {
    z <- c(x - stats::median(x), y - stats::median(y))
    zcomb <- utils::combn(z, 2)
    return(stats::median(abs(zcomb[1, ] - zcomb[2, ])))

  } else  if (type == "S3") {
    return(2 * stats::median(c(abs(x - stats::median(x)), abs(y - stats::median(y)))))

  } else if (type == "S4") {
    return(stats::median(abs(x - stats::median(x)) + stats::median(abs(y - stats::median(y)))))
  }
}

