## ----------------------------------------------------------------------------
## Location estimators for univariate samples and estimators for the location
## difference between two samples
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
#' ## Generate random sample
#' set.seed(108)
#' x <- rnorm(10)
#'
#' ## Compute 20% trimmed mean
#' trim_mean(x, gamma = 0.2)
#'
#' @export

trim_mean <- function(x, gamma = 0.2, na.rm = FALSE) {

  ## Check input arguments ----
  stopifnot(
    "'x' is missing." = !missing(x),
    "'x' must not be NULL." = !is.null(x),
    "'gamma' must not be NULL." = !is.null(gamma),
    "'na.rm' must not be NULL." = !is.null(na.rm),
    "'gamma' must not be NA." = !is.na(gamma),
    "'na.rm' must not be NA." = !is.na(na.rm),
    "'x' has to a numeric vector." = is.numeric(x),
    "'gamma' has to be a numeric value." = is.numeric(gamma),
    "'na.rm' has to be a logical value." = is.logical(na.rm),
    "'gamma' has to be a single value, not a vector of length >= 1." = identical(length(gamma), 1L),
    "'na.rm' has to be a single value, not a vector of length >= 1." = identical(length(na.rm), 1L)
  )

  if ((gamma < 0) || (gamma > 0.5)) {
    stop("'gamma' has to be a numeric value in [0, 0.5].")
  }

  ## Remove missing values in 'x' ----
  if (!na.rm & any(is.na(x))) {
    return(NA_real_)
  } else if (all(is.na(x))) {
    return(NA_real_)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Calculate trimmed mean ----
  return(mean(x, trim = gamma))
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
#' ## Generate random samples
#' set.seed(108)
#' x <- rnorm(10)
#'
#' ## Compute 20% winsorized mean
#' win_mean(x, gamma = 0.2)
#'
#' @export

win_mean <- function(x, gamma = 0.2, na.rm = FALSE) {

  ## Check input arguments ----
  stopifnot(
    "'x' is missing." = !missing(x),
    "'x' must not be NULL." = !is.null(x),
    "'gamma' must not be NULL." = !is.null(gamma),
    "'na.rm' must not be NULL." = !is.null(na.rm),
    "'gamma' must not be NA." = !is.na(gamma),
    "'na.rm' must not be NA." = !is.na(na.rm),
    "'x' has to a numeric vector." = is.numeric(x),
    "'gamma' has to be a numeric value." = is.numeric(gamma),
    "'na.rm' has to be a logical value." = is.logical(na.rm),
    "'gamma' has to be a single value, not a vector of length >= 1." = identical(length(gamma), 1L),
    "'na.rm' has to be a single value, not a vector of length >= 1." = identical(length(na.rm), 1L)
  )

  if ((gamma < 0) || (gamma > 0.5)) {
    stop("'gamma' has to be a numeric value in [0, 0.5].")
  }

  ## Remove missing values in 'x' ----
  if (!na.rm & any(is.na(x))) {
    return(NA_real_)
  } else if (all(is.na(x))) {
    return(NA_real_)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Calculate winsorized mean ----

  # For 'gamma == 0', the winsorized mean is identical to the sample mean
  if (identical(gamma, 0)) {
    return(mean(x))
  }

  # Number of trimmed observations
  n <- length(x)
  r <- floor(gamma * n)

  # Replace first and last r observations
  x.sort <- sort(x)
  x.sort[1:r] <- x.sort[r + 1]
  x.sort[(n - r + 1):n] <- x.sort[n - r]

  # Winsorized mean
  return(mean(x.sort))
}

#' @title One-sample Hodges-Lehmann estimator
#'
#' @description \code{hodges_lehmann} calculates the one-sample Hodges-Lehmann estimator
#' of a sample.
#'
#' @template x
#' @template na_rm
#'
#' @details The one-sample Hodges-Lehmann estimator for a sample of size \code{n}
#'          is defined as
#'
#' \deqn{med(\frac{X_i + X_j}{2},  1 \le i < j \le m).}
#'
#' @return
#' The one-sample Hodges-Lehmann estimator.
#'
#' @references
#' \insertRef{HodLeh63esti}{robTests}
#'
#' @examples
#' ## Generate random sample
#' set.seed(108)
#' x <- rnorm(10)
#'
#' ## Compute one-sample Hodges-Lehmann estimator
#' hodges_lehmann(x)
#'
#' @export

hodges_lehmann <- function(x, na.rm = FALSE) {

  ## Check input arguments ----
  stopifnot(
    "'x' is missing." = !missing(x),
    "'x' must not be NULL." = !is.null(x),
    "'na.rm' must not be NULL." = !is.null(na.rm),
    "'na.rm' must not be NA." = !is.na(na.rm),
    "'x' has to a numeric vector." = is.numeric(x),
    "'na.rm' has to be a logical value." = is.logical(na.rm),
    "'x' needs at least 2 values." = (length(x) > 1),
    "'na.rm' has to be a single value, not a vector of length >= 1." = identical(length(na.rm), 1L)
  )

  ## Remove missing values in 'x' ----
  if (!na.rm & any(is.na(x))) {
    return(NA_real_)
  } else if (all(is.na(x))) {
    return(NA_real_)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Calculate one-sample Hodges-Lehmann estimate ----

  # Compute pairwise means
  x.grid <- cbind(rep(1:length(x), each = length(x)), 1:length(x))
  x.diffs <- x.grid[x.grid[, 1] < x.grid[, 2], , drop = FALSE]
  mean.pairwise.sums <- (x[x.diffs[, 1]] + x[x.diffs[, 2]])/2

  ## Hodges-Lehmann estimate
  return(stats::median(mean.pairwise.sums))
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
#' @details The two-sample Hodges-Lehmann estimator for two samples \code{x}
#'          and \code{y} of sizes \code{m} and \code{n} is defined as
#'
#' \deqn{med(|x_i - y_j|, 1 \le i \le m, 1 \le j \le n).}
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
#' ## Generate random samples
#' set.seed(108)
#' x <- rnorm(10); y <- rnorm(10)
#'
#' ## Compute two-sample Hodges-Lehmann estimator
#' hodges_lehmann_2sample(x, y)
#'
#' @export

hodges_lehmann_2sample <- function(x, y, na.rm = FALSE) {

  ## Check input arguments ----
  stopifnot(
    "'x' is missing." = !missing(x),
    "'y' is missing." = !missing(y),
    "'x' must not be NULL." = !is.null(x),
    "'y' must not be NULL." = !is.null(y),
    "'na.rm' must not be NULL." = !is.null(na.rm),
    "'na.rm' must not be NA." = !is.na(na.rm),
    "'x' has to a numeric vector." = is.numeric(x),
    "'y' has to a numeric vector." = is.numeric(y),
    "'na.rm' has to be a logical value." = is.logical(na.rm),
    "'na.rm' has to be a single value, not a vector of length >= 1." = identical(length(na.rm), 1L)
  )

  # Remove missing values in 'x' and 'y' ----
  if (!na.rm & (any(is.na(x)) || any(is.na(y)))) {
    return(NA_real_)
  } else if (all(is.na(x)) || all(is.na(y))) {
    return(NA_real_)
  } else if (na.rm & (any(is.na(x)) || any(is.na(y)))) {
    x <- as.vector(stats::na.omit(x))
    y <- as.vector(stats::na.omit(y))
  }

  ## Calculate two-sample Hodges-Lehmann estimate ----

  # Compute pairwise differences between the two samples
  diff <- expand.grid(x, y)
  res <- diff[, 1] - diff[, 2]

  # Two-sample Hodges-Lehmann estimate
  return(stats::median(res))
}

#' @title M-estimator of location
#'
#' @description \code{m_est} calculates an M-estimate of location and its variance
#' for different psi functions.
#'
#' @template x
#' @template psi
#' @template k_mest
#' @template tol
#' @template max_it
#' @template na_rm
#'
#' @details
#' To compute the M-estimate, the iterative algorithm described in \insertCite{MarMarYoh06robu;textual}{robTests} is used.
#' The variance is estimated as in \insertCite{Hub81robu;textual}{robTests}.
#'
#' If \code{max.it} contains decimal places, it is truncated.
#'
#'
#' @return A list containing the components:
#'         \item{est}{estimated mean.}
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
#'
#' ## Generate random sample
#' set.seed(108)
#' x <- rnorm(10)
#'
#' ## Computer Huber's M-estimate
#' m_est(x, psi = "huber")
#'
#' @export

m_est <- function(x, psi, k = robustbase::.Mpsi.tuning.default(psi), tol = 1e-6, max.it = 15, na.rm = FALSE) {

  ## Check input arguments ----
  stopifnot(
    "'x' is missing." = !missing(x),
    "'psi' is missing." = !missing(psi),
    "'x' must not be NULL." = !is.null(x),
    "'psi' must not be NULL." = !is.null(psi),
    "'k' must not be NULL." = !is.null(k),
    "'tol' must not be NULL." = !is.null(tol),
    "'max.it' must not be NULL." = !is.null(max.it),
    "'na.rm' must not be NULL." = !is.null(na.rm),
    "'psi' must not be NA." = !is.na(psi),
    "'k' must not be NA." = !is.na(k),
    "'tol' must not be NA." = !is.na(tol),
    "'max.it' must not be NA." = !is.na(max.it),
    "'na.rm' must not be NA." = !is.na(na.rm),
    "'x' has to a numeric vector." = is.numeric(x),
    "'psi' has to a character value." = is.character(psi),
    "'k' has to be a numeric value." = is.numeric(k),
    "'tol' has to be a numeric value." = is.numeric(tol),
    "'max.it' has to be a numeric value." = is.numeric(max.it),
    "'na.rm' has to be a logical value." = is.logical(na.rm),
    "'psi' has to be a single value, not a vector of length >= 1." = identical(length(psi), 1L),
    "'psi' must be one of 'huber', 'hampel', or 'bisquare'." = any(c("huber", "hampel", "bisquare") %in% psi),
    # "'k' has to be a single value, not a vector of length >= 1." = identical(length(k), 1L),
    "'tol' has to be a single value, not a vector of length >= 1." = identical(length(tol), 1L),
    "'max.it' has to be a single value, not a vector of length >= 1." = identical(length(max.it), 1L),
    "'na.rm' has to be a single value, not a vector of length >= 1." = identical(length(na.rm), 1L),
    "'tol' has to be a positive value." = tol > 0,
    "'max.it' has to be a positive value." = max.it > 0
  )

  ## Remove missing values in 'x' ----
  if (!na.rm & any(is.na(x))) {
    return(NA_real_)
  } else if (all(is.na(x))) {
    return(NA_real_)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Remove decimal places from 'max.it' ----
  max.it <- trunc(max.it)

  ## Calculate M-estimate ----

  # Initial estimators
  est.old <- stats::median(x)
  S <- stats::mad(x)
  n.it <- 1

  # Iterative algorithm for computing the M-estimator, see. e.g. Maronna et al. (2006, p. 39)
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
