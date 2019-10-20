## ----------------------------------------------------------------------------
## Calculation of the different test statistics
## ----------------------------------------------------------------------------

#' @title Test statistic for the two-sample Yuen t-test
#'
#' @description
#' \code{trimmed_t} calculates the test statistic of the two-sample Yuen test.
#'
#' @template x
#' @template y
#' @template gamma_trimed_test
#' @template delta
#' @template na_rm
#'
#' @return
#' A list containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{estimates}{the trimmed means for both samples.}
#' \item{df}{the degrees of freedom for the test statistic.}
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertRef{YueDix73appr}{robTests}
#'
#' \insertRef{Yue74trim}{robTests}
#'
#' @export

trimmed_t <- function(x, y, gamma = 0.2, delta = 0, na.rm = FALSE) {

  ## Trimmed means
  x.trim <- trim_mean(x, gamma = gamma, na.rm = na.rm)
  y.trim <- trim_mean(y - delta, gamma = gamma, na.rm = na.rm)
  estimates <- c(x.trim, y.trim)

  ## Scale estimator
  var.x <- win_var(x, gamma = gamma, na.rm = na.rm)
  var.y <- win_var(y, gamma = gamma, na.rm = na.rm)

  h.x <- var.x$h
  h.y <- var.y$h
  var.x <- var.x$var
  var.y <- var.y$var

  m <- length(x)
  n <- length(y)
  pool.var <- ((m - 1) * var.x + (n - 1) * var.y)/(h.x + h.y - 2)

  ## Degrees of freedom
  df <- h.x + h.y - 2

  ## Test statistic
  statistic <- (x.trim - y.trim) / sqrt(pool.var * (1/h.x + 1/h.y))

  res <- list(statistic = statistic, estimates = estimates, df = df)

  return(res)
}

#' @title Test statistic for the asymmetrically trimmed test
#'
#' @description
#' \code{asym_trimmed_t} calculates the test statistic of the asymmetrically trimmed test with a given skewness selector statistic
#'
#' @template x
#' @template y
#' @template type_skewness
#' @template delta
#' @template na_rm
#'
#' @return
#' A list containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{estimates}{the asymmetrically trimmed means for both samples.}
#' \item{df}{the degrees of freedom for the test statistic.}
#'
#' @references
#' \insertRef{YueDix73appr}{robTests}
#'
#' \insertRef{Yue74trim}{robTests}
#'
#' @export

asym_trimmed_t <- function(x, y, type, delta = 0, na.rm = FALSE) {

  ## Trimmed means
  x.trim <- asym_trimmed_mean(x, type = type, na.rm = na.rm)
  y.trim <- asym_trimmed_mean(y - delta, type = type, na.rm = na.rm)
  estimates <- c(x.trim, y.trim)

  ## Scale estimator
  var.x <- asym_win_var(x, type = type, na.rm = na.rm)
  var.y <- asym_win_var(y, type = type, na.rm = na.rm)

  h.x <- var.x$h
  h.y <- var.y$h
  var.x <- var.x$var
  var.y <- var.y$var

  m <- length(x)
  n <- length(y)
  pool.var <- ((h.x - 1) * var.x + (h.y - 1) * var.y)/(h.x + h.y - 2)

  ## Degrees of freedom
  df <- h.x + h.y - 2

  ## Test statistic
  statistic <- (x.trim - y.trim) / sqrt(pool.var * (1/h.x + 1/h.y))

  res <- list(statistic = statistic, estimates = estimates, df = df)

  return(res)
}


#' @title Robust permutation statistics based on medians
#'
#' @description \code{rob_perm_statistic()} calculates test statistics for robust permutation tests based on medians.
#'
#' @template x
#' @template y
#' @template type_rob_perm
#' @template na_rm
#'
#' @return A list containing the following components:
#'         \item{statistic}{the selected test statistic}
#'         \item{abs.statistic}{the absolute value of the selected test statistic}
#'         \item{estimates}{estimate of location for each sample if available}
#'
#' @details The test statistics returned by \code{rob_perm_statistic} are of the form
#'          \deqn{D_i/S_j} where the D_i, i = 1,...,3, are different estimators of location and the S_j, j = 1,...,4 are
#'          estimates for the mutual sample scale. See Fried and Dehling (2011) for details.
#'
#' @examples
#' x <- rnorm(20); y <- rnorm(20)
#' rob_perm_statistic(x, y, type = "D2S1")
#'
#' @references
#' \insertRef{FriDeh11robu}{robTests}
#'
#' @export

rob_perm_statistic <- function(x, y,
                          type = c("D1S1", "D1S2", "D2S1", "D2S2", "D3S3", "D3S4"),
                          na.rm = FALSE) {
  type <- match.arg(type)

  if(!(type %in% c("D1S1", "D1S2", "D2S1", "D2S2", "D3S3", "D3S4"))) {
    stop("type needs to be one of 'D1S1', 'D1S2', 'D2S1', 'D2S2', 'D3S3', 'D3S4'")
  }

  switch(type,
         D1S1 = {
            est.x <- hodges_lehmann(x, na.rm = na.rm)
            est.y <- hodges_lehmann(y, na.rm = na.rm)
            loc <- est.x - est.y
            sd <- rob_var(x, y, na.rm = na.rm, type = "S1")
            res <- loc/sd
         },
         D1S2 = {
            est.x <- hodges_lehmann(x, na.rm = na.rm)
            est.y <- hodges_lehmann(y, na.rm = na.rm)
            loc <- est.x - est.y
            sd <- rob_var(x, y, na.rm = na.rm, type = "S2")
            res <- loc/sd
         },
         D2S1 = {
            est.x <- est.y <- NULL
            loc <- hodges_lehmann_2sample(x, y, na.rm = na.rm)
            sd <- rob_var(x, y, na.rm = na.rm, type = "S1")
            res <- loc/sd
         },
         D2S2 = {
            est.x <- est.y <- NULL
            loc <- hodges_lehmann_2sample(x, y, na.rm = na.rm)
            sd <- rob_var(x, y, na.rm = na.rm, type = "S2")
            res <- loc/sd
         },
         D3S3 = {
            est.x <- stats::median(x, na.rm = na.rm)
            est.y <- stats::median(y, na.rm = na.rm)
            loc <- est.x - est.y
            sd <- rob_var(x, y, na.rm = na.rm, type = "S3")
            res <- loc/sd
        },
         D3S4 = {
            est.x <- stats::median(x, na.rm = na.rm)
            est.y <- stats::median(y, na.rm = na.rm)
            loc <- est.x - est.y
            sd <- rob_var(x, y, na.rm = na.rm, type = "S4")
            res <- loc/sd
      })

  return(list(statistic = res, abs.statistic = abs(res), estimates = c(est.x, est.y)))
}


#' @title Calculate the M-test-statistic
#'
#' @description \code{m_test_statistic} calculates the test statistics for randomization tests based on M-estimators
#'
#' @template x
#' @template y
#' @template k_mest
#'
#' @return A list containing
#'         \item{statistic}{the standardized test statistic, and}
#'         \item{estimates}{the M-estimates of location for both \code{x} and \code{y}}
#'
#' @export

m_test_statistic <- function(x, y, psi, k = .Mpsi.tuning.default(psi)) {
  ## Sample sizes
  m <- length(x)
  n <- length(y)

  ## M-estimators and corresponding variances for both samples
  m.x <- m_est(x, psi = psi, k = k, max.it = 1)
  m.y <- m_est(y, psi = psi, k = k, max.it = 1)

  est.x <- m.x$est
  est.y <- m.y$est

  var.x <- m * m.x$var
  var.y <- n * m.y$var

  ## Test statistic
  return(list(statistic = (est.x - est.y) / sqrt(((m - 1) * var.x + (n - 1) * var.y) / (n + m - 2) * (1/m + 1/n)),
              estimates = c(est.x, est.y)))
}


#' @title Simultaneous Huber-M-estimates of scale and location
#'
#' @description Calculates M-estimates of location and the joined scale of two samples
#'
#' @template x
#' @template y
#' @template k_mest
#'
#' @return Named list containing the following objects
#'         \item{mu.x}{Location estimate of x}
#'         \item{mu.y}{Location estimate of y}
#'         \item{s}{Scale estimate for the joined sample}
#' @import robustbase
#' @export

huber_2sample <- function(x, y, k) {
  m <- length(x)
  n <- length(y)

  N <- m + n

  beta <- 2 * k^2 * (1 - stats::pnorm(k)) + 2 * stats::pnorm(k) - 1 - sqrt(2/pi) * k * exp(-1/2 * k^2)

  s.old <- stats::mad(x) + stats::mad(y)
  #2 * stats::median(c(abs(x - stats::median(x)), abs(y - stats::median(y))))
  mux.old <- stats::median(x)
  muy.old <- stats::median(y)

  #repeat {
  z.x <- (x - mux.old)/s.old
  z.y <- (y - muy.old)/s.old

  s.new <- sqrt(1/((N - 1) * beta) * (sum(robustbase::Mpsi(z.x, psi = "huber", cc = k)^2) + sum(robustbase::Mpsi(z.y, psi = "huber", cc = k)^2)) * s.old^2)

  mux.new <- mux.old + (1/m * sum(robustbase::Mpsi(z.x, psi = "huber", cc = k)) * s.old)/(1/m * sum(robustbase::Mpsi(z.x, psi = "huber", cc = k, deriv = 1)))

  muy.new <- muy.old + (1/n * sum(robustbase::Mpsi(z.y, psi = "huber", cc = k)) * s.old)/(1/n * sum(robustbase::Mpsi(z.y, psi = "huber", cc = k, deriv = 1)))


  #if (abs(mux.new - mux.old) < 1e-6 & abs(muy.new - muy.old) < 1e-6 & abs(s.new - s.old) < 1e-6) {
  #  break
  #}

  s.old <- s.new
  mux.old <- mux.new
  muy.old <- muy.new
  #}

  return(list(mu.x = mux.new, mu.y = muy.new, s = s.new))
}

