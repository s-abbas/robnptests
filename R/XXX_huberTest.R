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


## ----------------------------------------------------------------------------
## Huber two-sample test
## ----------------------------------------------------------------------------

#' @title Two-sample location tests based on Huber's M-estimator
#'
#' @description Performs a 2-sample location test based on Huber's M-estimator
#'
#' @template x
#' @template y
#' @template delta
#' @template k_huber
#' @template alternative
#' @template var_test
#' @template na_rm
#'
#' @return
#' A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom for the test statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{estimate}{the Huber M-estimates of \code{x} and \code{y}.}
#' \item{null.value}{the specified hypothesized value of the mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating how the p-value was computed.}
#' \item{data.name}{a character string giving the names of the data.}
#'
#' @details
#' The test is introduced in the context of the hybrid tests in Weichert & Hothorn (2002).
#'
#' @examples
#' x <- rnorm(20); y <- rnorm(20)
#' huber_test(x, y)
#'
#' @seealso
#' \code{\link[robTests]{min_tc_test}}
#' \code{\link[robTests]{min_t_test}}
#'
#' @importFrom Rdpack reprompt
#'
#' @references \insertRef{WeiHot02robu}{robTests}
#'
#' @export

huber_test <- function(x, y, delta = ifelse(var.test, 1, 0), k = 1.8,
                       alternative = c("two.sided", "greater", "less"),
                       var.test = FALSE, na.rm=FALSE) {

  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  ## NA handling
  if (!na.rm & (any(is.na(x)) | any(is.na(y)))) {
    return(NA)
  } else if (na.rm & (any(is.na(x)) | any(is.na(y)))) {
    x <- as.numeric(stats::na.omit(x))
    y <- as.numeric(stats::na.omit(y))
  }

  ## If necessary: Transformation to test for difference in scale
  if (var.test) {
    x <- log(x^2)
    y <- log(y^2)
    delta <- log(delta^2)
  }

  alternative <- match.arg(alternative)

  m <- length(x)
  n <- length(y)

  hub <- huber_2sample(x, y - delta, k = k)

  m.x <- hub$mu.x
  m.y <- hub$mu.y

  stand <- hub$s

  #statistic <- (m.x$mu - m.y$mu)/sqrt((1/n.x + 1/n.y) * ((n.x - 1) * m.x$s^2 + (n.y - 1) * m.y$s^2)/(n.x + n.y - 2))
  statistic <- (m.x - m.y)/(sqrt(1/m + 1/n) * stand)

  if (delta != 0) {
    estimates <- c(m.x, huber_2sample(x, y, k = k)$mu.y)
  } else estimates <- c(m.x, m.y)

  p.value <- switch (alternative,
                     two.sided = 2 * stats::pt(abs(statistic), df = m + n - 2, lower.tail = FALSE),
                     greater = stats::pt(statistic, df = m + n - 2, lower.tail = FALSE),
                     less = stats::pt(statistic, df = m + n - 2, lower.tail = TRUE)
  )


  if (var.test) {
    names(estimates) <- c("Huber-M of log(x^2)", "Huber-M of log(y^2)")
    names(delta) <- "ratio of variances"
    delta <- exp(delta)
  } else {
    names(estimates) <- c("Huber-M of x", "Huber-M of y")
    names(delta) <- "location shift"
  }

  names(statistic) <- ifelse(var.test, "S", "D")
  method <- "Huber-Test"
  parameter <- m + n - 2
  names(parameter) <- "df"

  res <- list(statistic = statistic, parameter = parameter, p.value = p.value,
              estimate = estimates, null.value = delta,
              alternative = alternative,
              method = method, data.name = dname)

  class(res) <- "htest"

  return(res)
}



