## ----------------------------------------------------------------------------
## Location estimators for univariate samples and
## estimators for the location difference between two samples
## ----------------------------------------------------------------------------

#' @title Trimmed mean
#'
#' @description
#' \code{trim_mean} calculates a trimmed mean of a sample.
#'
#' @template x
#' @template gamma_trimmed_mean
#' @template na_rm
#'
#' @details
#' This is a wrapper function for the function \code{\link[base]{mean}}.
#'
#' @return
#' The trimmed mean.
#'
#' @examples
#' set.seed(108)
#' x <- rnorm(10)
#' trim_mean(x, gamma = 0.2)
#'
#' @export

trim_mean <- function(x, gamma = 0.2, na.rm = FALSE) {
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

  ## Calculate trimmed mean
  return(mean(x, trim = gamma))
}

#' @title Asymmetrically trimmed mean
#'
#' @description
#' \code{asym_trimmed_mean} calculates an asymmetrically trimmed mean of a sample.
#' In contrast to an ordinary trimmed mean, the numbers of observations removed
#' from the lower and the upper end of the sample do not need to be equal.
#'
#' @template x
#' @template type_skewness
#' @template na_rm
#'
#' @details
#' The number of observations trimmed from the lower and the upper end of the sample
#' depends on the skewness.
#'
#' Three skewness-selector statistics suggested by Reed and Stark (1996) have been
#' implemented. The argument \code{type} specifies which one is used.
#'
#' \describe{
#' \item{"Q2"}{This skewness-selector statistic compares the difference between the upper and the lower 5\% of
#'             the observations to the difference between the upper and lower 50\%.}
#' \item{"SK2"}{This skewness-selector statistic compares the difference between the sample minimum and
#'              median to the difference between the median and the maximum.}
#' \item{"SK5"}{This skewness-selector statistic compares the difference between the sample minimum and
#'              mean to the difference between the mean and the maximum.}
#' }
#'
#' @return
#' The asymmetrically trimmed mean.
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{ReeSta96hing}{robTests}
#'
#' @examples
#' set.seed(108)
#' x <- rnorm(10)
#' asym_trimmed_mean(x, type = "Q2")
#' asym_trimmed_mean(x, type = "SK2")
#' asym_trimmed_mean(x, type = "SK5")
#'
#' @export

asym_trimmed_mean <- function(x, type = c("Q2", "SK2", "SK5"), na.rm = FALSE) {
  ## NA handling
  if (!na.rm & any(is.na(x))) {
    return(NA)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  type <- match.arg(type)

  ## Sample size and sorted sample
  m <- length(x)
  x.sort <- sort(x)

  ## Determine lower trimming proportion
  if (type == "Q2") {

    ## Do not trim, if sample size is too small
    if (floor(0.05 * m) == 0) {
      return(mean(x))
    }

    gamma <- 0.1

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

  ## Calculate number of values to be trimmed from both ends
  r.lower <- floor(m * gamma.lower) + 1
  r.upper <- m - floor(m * gamma.upper)

  ## Asymmetrically trimmed mean
  res <- mean(x.sort[(r.lower + 1):r.upper])

  return(res)
}

#' @title Winsorized mean
#'
#' @description \code{win_mean} calculates the winsorized mean of a sample.
#'
#' @template x
#' @template gamma_winsorized_mean
#' @template na_rm
#'
#' @return
#' The winsorized mean.
#'
#' @examples
#' set.seed(108)
#' x <- rnorm(10)
#' win_mean(x, gamma = 0.2)
#'
#' @export

win_mean <- function(x, gamma = 0.2, na.rm = FALSE) {
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

  n <- length(x)

  ## Number of trimmed observations
  r <- floor(gamma * n)

  ## Replace first and last r observations
  x.sort <- sort(x)
  x.sort[1:r] <- x.sort[r + 1]
  x.sort[(n - r + 1):n] <- x.sort[n - r]

  ## Winsorized mean
  res <- mean(x.sort)

  return(res)
}

#' @title One-sample Hodges-Lehmann estimator
#'
#' @description \code{hodges_lehmann} calculates the one-sample Hodges-Lehmann estimator
#' of a sample.
#'
#' @template x
#' @template na_rm
#'
#' @return
#' The one-sample Hodges-Lehmann estimator.
#'
#' @references
#' \insertRef{HodLeh63esti}{robTests}
#'
#' @examples
#' set.seed(108)
#' x <- rnorm(10)
#' hodges_lehmann(x)
#'
#' @export

hodges_lehmann <- function(x, na.rm = FALSE) {

  ## NA handling
  if (!na.rm & any(is.na(x))) {
    return(NA)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Pairwise means
  x.grid <- cbind(rep(1:length(x), each = length(x)), 1:length(x))
  x.diffs <- x.grid[x.grid[, 1] < x.grid[, 2], , drop = FALSE]
  mean.pairwise.sums <- (x[x.diffs[, 1]] + x[x.diffs[, 2]])/2

  ## Hodges-Lehmann estimate
  res <- stats::median(mean.pairwise.sums)

  return(res)
}


#' @title Two-sample Hodges-Lehmann estimator
#'
#' @description \code{hodges_lehmann_2sample} calculates the two-sample Hodges-Lehmann
#' estimator for the location difference of two samples x and y.
#'
#' @template x
#' @template y
#' @template na_rm
#'
#' @return
#' The two-sample Hodges-Lehmann estimator.
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{HodLeh63esti}{robTests}
#'
#' @examples
#' set.seed(108)
#' x <- rnorm(10); y <- rnorm(10)
#' hodges_lehmann_2sample(x, y)
#'
#' @export

hodges_lehmann_2sample <- function(x, y, na.rm = FALSE) {

  ## NA handling
  if (!na.rm & (any(is.na(x)) || any(is.na(y)))) {
    return(NA)
  } else if (na.rm & (any(is.na(x)) || any(is.na(y)))) {
    x <- as.vector(stats::na.omit(x))
    y <- as.vector(stats::na.omit(y))
  }

  diff <- expand.grid(x, y)
  res <- diff[, 1] - diff[, 2]

  return(stats::median(res))
}

#' @title M-estimator of location
#'
#' @description \code{m_est} calculates an M-estimate of location and its variance using different tuning functions.
#'
#' @template x
#' @template psi
#' @template k_mest
#' @template tol
#' @template max_it
#' @template na_rm
#'
#' @details
#' To compute the M-estimate, the iterative algorithm described in Maronna et al. (2006, p. 39) is used.
#' The variance is estimated according to Huber (1981, p. 150).
#'
#'
#' @return A list containing the items:
#'         \item{est}{estimated mean, and}
#'         \item{var}{estimated variance.}
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{MarMarYoh06robu}{robTests}
#'
#' \insertRef{Hub81robu}{robTests}
#'
#'
#' @examples
#' set.seed(108)
#' x <- rnorm(10)
#'
#' m_est(x, psi = "huber")
#'
#' @export

m_est <- function(x, psi, k = robustbase::.Mpsi.tuning.default(psi), tol = 1e-6, max.it = 15, na.rm = TRUE) {

  ## NA handling
  if (!na.rm & any(is.na(x))) {
    return(NA)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Initial estimators
  est.old <- stats::median(x)
  S <- stats::mad(x)
  n.it <- 1

  ## Iterative algorithm for computing the M-estimator, see. e.g. Maronna et al. (2006, p. 39)
  repeat {
    w <- robustbase::Mwgt((x - est.old)/S, psi = psi, cc = k)
    est.new <- sum(w * x)/sum(w)

    if (abs(est.new - est.old) < tol | n.it >= max.it) {
      break
    }

    est.old <- est.new
    n.it <- n.it + 1
  }

  est <- est.new

  ## Variance estimation, see e.g. Huber (1981, p. 150)
  n <- length(x)
  z <- (x - est)/S

  var <- 1/(n * (n - 1)) * S^2 * sum(robustbase::Mpsi(z, psi = psi, cc = k)^2)/(1/n * sum(robustbase::Mpsi(z, psi = psi, cc = k, deriv = 1)))^2

  return(list(est = est, var = var))
}
