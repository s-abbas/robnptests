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
  checkmate::assert_numeric(x, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_number(gamma, na.ok = FALSE, lower = 0, upper = 0.5, finite = TRUE, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)

  ## Remove missing values in 'x' ----
  if (!na.rm & any(is.na(x))) {
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
  checkmate::assert_numeric(x, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_number(gamma, na.ok = FALSE, lower = 0, upper = 0.5, finite = TRUE, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)

  ## Remove missing values in 'x' ----
  if (!na.rm & any(is.na(x))) {
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
  checkmate::assert_numeric(x, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)

  ## Remove missing values in 'x' ----
  if (!na.rm & any(is.na(x))) {
    return(NA_real_)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

  ## Calculate one-sample Hodges-Lehmann estimate ----

  # Compute pairwise means
  x.grid <- cbind(rep(seq_along(x), each = length(x)), seq_along(x))
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
  checkmate::assert_numeric(x, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(y, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)

  # Remove missing values in 'x' and 'y' ----
  if (!na.rm & (any(is.na(x)) || any(is.na(y)))) {
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
#' To compute the M-estimate, the iterative algorithm described in
#' \insertCite{MarMarYoh06robu;textual}{robTests} is used.
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
  checkmate::assert_numeric(x, finite = TRUE, all.missing = FALSE, null.ok = FALSE)
  checkmate::assert_choice(psi, choices = c("huber", "hampel", "bisquare"), null.ok = FALSE)
  checkmate::assert_numeric(k, lower = 0, len = ifelse(psi == "hampel", 3, 1), finite = TRUE, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_number(tol, na.ok = FALSE, lower = 1e-6, finite = TRUE, null.ok = FALSE)
  checkmate::assert_count(max.it, na.ok = FALSE, positive = TRUE, null.ok = FALSE)
  checkmate::assert_flag(na.rm, na.ok = FALSE, null.ok = FALSE)

  ## Remove missing values in 'x' ----
  if (!na.rm & any(is.na(x))) {
    return(NA_real_)
  } else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }

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
